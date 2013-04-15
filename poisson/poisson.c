#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <mpi.h>
#include <omp.h>

typedef struct{
    int local_m;
    int m;
    double** data;
    int comm_size;
    int comm_rank;
    int* sizes;
    int* displ;
    int* count;
} Matrix;

/* Function prototypes */
double *createVector(int n);
double **createMatrix(int m, int n);
void freeMatrix(Matrix x);
void fst_(double *v, int *n, double *w, int *nn);
void fstinv_(double *v, int *n, double *w, int *nn);
void transposeMPI(Matrix ut, Matrix u);
double evalFunc(int i, int j, double h, int displ, double pi);
void printToFile(double**u, int n, int m, int rank);

int main (int argc, char** argv){

    if (argc < 2) {
        printf("Usage: %s <N> [L]\n",argv[0]);
        return 1;
    }

    const double pi = 4.0*atan(1.0);
    double *lambda, **z, *globalDispl;
    Matrix u;
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int n = atoi(argv[1]);
    int m = n-1;
    int local_m = m/size;
    int nn = 4*n;
    double len = 1.0;
    if (argc > 2)
        len = atof(argv[2]);
    double h = len/(double)n;
    double hh = h*h;

    int threads = omp_get_max_threads();

    lambda = createVector(m);
    globalDispl = createVector(size+1); //+1 for sleeker code
    z = createMatrix(threads,nn);

    u.count = calloc(size, sizeof(int));
    u.sizes = calloc(size, sizeof(int));
    u.displ = calloc(size, sizeof(int));
    u.m = m;
    u.comm_size = size;
    u.comm_rank = rank;


    //Calculate number of rows each processor should own.
    for (int i = 0; i < size; i++){
        u.sizes[i] = local_m;
        if (i < m%local_m)
            u.sizes[i] += 1;
        globalDispl[i+1] = globalDispl[i] + u.sizes[i];
    }
    if ((local_m*size < m) && (rank < m%local_m))
        local_m++;

    //Calculate how many elts. to send/recv and processor displacement.
    int sum = 0;
    for (int i = 0; i < size; i++){
        u.count[i] = u.sizes[rank]*u.sizes[i];
        u.displ[i] = sum;
        sum += u.count[i];
    }

    u.data = createMatrix(local_m+1, m);

    //Calculate eigenvalues.
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < m; i++){
        lambda[i] = 2.0*(1.0-cos((i+1)*pi*h));
    }

    //Initialize G
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < local_m; i++){
        for (int j = 0; j < m; j++){
            u.data[i][j] = hh*evalFunc(i, j, h, globalDispl[rank], pi);
        }
    }
    double startTime = MPI_Wtime();

    //Step 1) G~t = QtGQ = S⁻¹((S(G))t)
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < local_m; i++) {
        fst_(u.data[i], &n, z[omp_get_thread_num()], &nn);
    }
    transposeMPI(u,u);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < local_m; i++){
        fstinv_(u.data[i], &n, z[omp_get_thread_num()], &nn);
    }

    //Step 2) ũ_ji = g~_ji / (l_j + l_i)
    #pragma omp parallel for schedule(static)
    for (int j=0; j < local_m; j++){
        for (int i=0; i < m; i++){
            u.data[j][i] /= (lambda[(int)globalDispl[rank]+j]+lambda[i]);
        }
    }

    //Step 3) U = QŨQ^T = S⁻¹((S(Ut))t)
    #pragma omp parallel for schedule(static)
    for (int i=0; i < local_m; i++){
        fst_(u.data[i], &n, z[omp_get_thread_num()], &nn);
    }
    transposeMPI(u, u);

    #pragma omp parallel for schedule(static)
    for (int i=0; i < local_m; i++){
        fstinv_(u.data[i], &n, z[omp_get_thread_num()], &nn);
    }

    double umax = 0.0, temp;
    for (int i = 0; i < local_m; i++){
        for (int j = 0; j < m; j++){
            temp = u.data[i][j]-evalFunc(i, j, h, globalDispl[rank], pi)/(5.*pi*pi);
            if (temp > umax){
                umax = temp;
            }
        }
    }
    double totTime = MPI_Wtime()-startTime;
    double t = totTime;
    double max;
    MPI_Reduce(&t, &totTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Allreduce(&umax, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (rank == 0){
	printf("Processors (MPI,threads):     %i*%i=%i\n",size,threads,size*threads);
        printf("Jobsize:                      %i/%i\n",local_m,m);
        printf("Average runtime:              %.8f s\n", totTime/size);
        printf("Max pointwise error:          %.14f \n\n", max);
    }

    //printToFile(u.data, local_m, m, rank);

    free(u.data);
    free(lambda);
    free(z);
    free(globalDispl);
    MPI_Finalize();

    return 0;
}

double *createVector (int n){
	double *a;
	int i;
	a = (double*)malloc(n*sizeof(double));
	for (i=0; i < n; i++) {
		a[i] = 0.0;
	}
	return a;
}

double **createMatrix(int n1, int n2){
	double **a;
	a    = (double **)malloc(n1   *sizeof(double *));
	a[0] = (double  *)malloc(n1*n2*sizeof(double));
	for (int i=1; i < n1; i++) {
		a[i] = a[i-1] + n2;
	}
	return (a);
}

double evalFunc(int i, int j, double h, int displ, double pi){
	double x = (double)(j+1)*h;
	double y = (double)(displ+i+1)*h;
	//printf("(%1.8f, %1.8f)\t", x,y);
	//if (x > 0.4 && y > 0.4){return 15.;}
	//if (x == 0.25 && y == 0.25){printf("(%1.8f, %1.8f) = 10\n", x,y);return -5.;}
	//if (x == 0.75 && y == 0.75){printf("(%1.8f, %1.8f) = -10\n", x,y);return 8.;}
	//return exp(1.-x-y)-1.;
	//return -(2.0-4*pi*pi*x*(x-1.0))*sin(2.0*pi*y);
	return 5.0*pi*pi*sin(pi*x)*sin(2.0*pi*y);
	//return 1.0;

}

void transposeMPI(Matrix ut, Matrix u){
	int len = ut.displ[ut.comm_size-1]+ut.count[ut.comm_size-1];
	double* sendbuff = (double*)malloc(len*sizeof(double));
	double* recvbuff = (double*)malloc(len*sizeof(double));
	int l = 0;
	for (int i = 0,count = 0; i < ut.comm_size; i++){
		for (int j = 0; j < ut.sizes[ut.comm_rank]; j++){
			for (int k = 0; k < ut.sizes[i]; k++){
				sendbuff[count++] = u.data[j][k+l];
			}
		}
		l += ut.sizes[i];
	}
	MPI_Alltoallv(sendbuff,u.count,u.displ,MPI_DOUBLE,recvbuff,ut.count,ut.displ,MPI_DOUBLE,MPI_COMM_WORLD);
	free(sendbuff);
	for (int i = 0,count = 0; i < ut.m; i++){
		for (int j = 0; j < ut.sizes[ut.comm_rank]; j++){
			ut.data[j][i] = recvbuff[count++];
		}
	}
	free(recvbuff);
}

void printToFile(double**u, int n, int m, int rank){
    printf("starting write\n");
	FILE *file;
	char filename[16];
	sprintf(filename, "out%d.dat", rank);
	file = fopen(filename, "w+");

	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
			fprintf(file, "%.10f ", u[i][j]);
		}
		fprintf(file, "\n");
	}
	fclose(file);
}

































