#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <mpi.h>
#include <omp.h>

typedef struct{
    int m;
    int local_m;
    double** data;
    int comm_size;
    int comm_rank;
    int* sizes;
    int* displ;
    int* count;
} Matrix;

/* Prototypes */
double *createVector(int n);
double **createMatrix(int m, int n);
void fst_(double *v, int *n, double *w, int *nn);
void fstinv_(double *v, int *n, double *w, int *nn);
void transposeMPI(Matrix ut, Matrix u);
double evalFunc(int i, int j, double h, int displ, double pi);

int main (int argc, char** argv){

    if (argc < 2) {
        printf("Usage: %s <N> [L]\n",argv[0]);
        return 1;
    }

    const double pi = 4.0*atan(1.0);
    double *lambda, **z, *globalDispl;
    Matrix u, ut;
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

    u.count = ut.count = calloc(size, sizeof(int));
    u.sizes = ut.sizes = calloc(size, sizeof(int));
    u.displ = ut.displ = calloc(size, sizeof(int));
    u.m = ut.m = m;
    u.comm_size = ut.comm_size = size;
    u.comm_rank = ut.comm_rank = rank;

    //Calculate number of rows each processor should have.
    for (int i = 0; i < size; i++){
        u.sizes[i] = ut.sizes[i] = local_m;
        if (i < m%local_m)
            u.sizes[i] = ut.sizes[i] += 1;
        globalDispl[i+1] = globalDispl[i] + u.sizes[i];
    }

    //Calculate how many elts. to send/recv and processor displ.
    int sum = 0;
    for (int i = 0; i < size; i++){
        u.count[i] = ut.count[i] = u.sizes[rank]*u.sizes[i];
        //printf("%i Rank: %i \t Count: %i\n",i,rank,u.count[i]);
        u.displ[i] = ut.displ[i] = sum;
        sum = sum + u.count[i]; //sum += u.count[i];
    }

    //Calculate local number of rows.
    if ((local_m*size < m) && (rank < m%local_m)) // <<-- m%local_m
        local_m++;

    u.data = createMatrix(local_m, m);
    ut.data = createMatrix(local_m, m);

    //Calculate eigenvalues.
    for (int i = 0; i < m; i++){
        lambda[i] = 2.0*(1.0-cos((i+1)*pi*h));
    }
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < local_m; i++){
        for (int j = 0; j < m; j++){
            u.data[i][j] = hh*evalFunc(i, j, h, globalDispl[rank], pi);
        }
    }
    double startTime = MPI_Wtime();
/*if (rank == 0){
printf("Init\n");
printf("data: %f\n",u.data[0][0]);
}*/
    //Step 1) G = QtGQ = S⁻¹((S(G))t)
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < local_m; i++) {
        fst_(u.data[i], &n, z[omp_get_thread_num()], &nn);
    }

    transposeMPI(ut,u);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < local_m; i++){
        fstinv_(ut.data[i], &n, z[omp_get_thread_num()], &nn);
    }
/*if (rank == 0){
printf("Step 1\n");
printf("data: %e\n",ut.data[0][0]);
}*/
    //Step 2) u_ij = g_ij / (l_i + l_j)
    #pragma omp parallel for schedule(static)
    for (int i=0; i < local_m; i++){
        for (int j=0; j < m; j++){
            ut.data[i][j] /= (lambda[(int)globalDispl[rank]+i]+lambda[j]);//lambda[(int)globalDispl[rank]+i]+lambda[j];
        }
    }
/*if (rank == 0){
printf("Step 2\n");
printf("data: %e\n",ut.data[3][3]);
}*/
    //Step 3) U = QUQ^T = S⁻¹((S(Ut))t)
    #pragma omp parallel for schedule(static)
    for (int i=0; i < local_m; i++){
        fst_(ut.data[i], &n, z[omp_get_thread_num()], &nn);
    }
    transposeMPI(u, ut);

    #pragma omp parallel for schedule(static)
    for (int i=0; i < local_m; i++){
        fstinv_(u.data[i], &n, z[omp_get_thread_num()], &nn);
    }

    double umax = 0.0, temp;
    for (int i = 0; i < local_m; i++){
        for (int j = 0; j < m; j++){
            temp = u.data[i][j]-evalFunc(i, j, h, globalDispl[rank], pi)/(5.*pi*pi);
            //printf("%i,%i data: %f\n",i,j,u.data[i][j]);
            if (temp > umax){
                umax = temp;
            }
        }
    }

    MPI_Allreduce(&umax, &umax, 0, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (rank == 0)
        printf ("Rank: %i Max pointwise error = %.14f \n", rank, umax);

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
	//memset(a[0],0,n1*n2*sizeof(double));
	return (a);
}

double evalFunc(int i, int j, double h, int displ, double pi){
    double x = (double)(j+1)*h;
    double y = (double)(displ+i+1)*h;
    //if (x > 0.6 && y > 0.6) return 10.;
    //return exp(1.-x-y)-1.;
    //return -(2.0-4*pi*pi*x*(x-1.0))*sin(2.0*pi*y);
	return 5.0*pi*pi*sin(pi*x)*sin(2.0*pi*y);
    //return 1.0;

}

void transposeMPI(Matrix ut, Matrix u){
	int i,j,k;
	int len = ut.displ[ut.comm_size-1]+ut.count[ut.comm_size-1];
	double* temp = (double*)malloc(len*sizeof(double));
	double* temp2 = (double*)malloc(len*sizeof(double));
	int l = 0;
	int count = 0;
	for (i = 0; i < ut.comm_size; i++){
		for (j = 0; j < ut.sizes[ut.comm_rank]; j++){
			for (k = 0; k < ut.sizes[i]; k++){
				temp[count++] = u.data[j][k+l];
			}
		}
		l += ut.sizes[i];
	}

	MPI_Alltoallv(temp,u.count,u.displ,MPI_DOUBLE,temp2,ut.count,ut.displ,MPI_DOUBLE,MPI_COMM_WORLD);
	count = 0;
	for (i = 0; i < ut.m; i++){
		for (j = 0; j < ut.sizes[ut.comm_rank]; j++){
			ut.data[j][i] = temp2[count++];
		}
	}
	free(temp);
	free(temp2);
}

































