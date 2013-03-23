#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <mpi.h>
#include <omp.h>

#ifndef HAVE_MKL
#define fst fst_
#define fstinv fstinv_
#endif


void fst(double *v, int *n, double *w, int *nn);
void fstinv(double *v, int *n, double *w, int *nn);
typedef struct {
    int m;
    int M;
    double** data;
    int comm_size;
    int comm_rank;
    int* sizes;
    int* displ;
    int* count;
} Matrix;
double findMax(Matrix u){
    double max = 0.0;
    for (int i=0; i < u.m; i++) {
        for (int j=0; j < u.M; j++) {
            if (u.data[j][i] > max) max = u.data[j][i];
        }
    }
  return max;
}
double **createMatrix (int n1, int n2){
	double **a;
	a    = (double **)malloc(n1*sizeof(double *));
	a[0] = (double *)malloc(n1*n2*sizeof(double));
	for (int i=1; i < n1; i++) {
		a[i] = a[i-1] + n2;
	}
	int n = n1*n2;
	memset(a[0],0,n*sizeof(double));
	return (a);
}
double *createVector (int n){
	double *a;
	a = (double*)malloc(n*sizeof(double));
	for (int i=0; i < n; i++) {
		a[i] = 0.0;
	}
	return a;
}
/*double *createEigenvalues(int m, double pi){
  double* diag = createVector(m);
  for (int i=0; i < m; i++)
    diag[i] = 2.0*(1.0-cos((i+1)*pi/(m+1)));
  return diag;
}*/
double funcEval(int i, int j, double h, int displ, double pi){
	double y = (double)(displ+i+1)*h;
	double x = (double)(j+1)*h;
    //if (x > 0.6 && y > 0.6) return 10.;
    //return exp(1.-x-y)-1.;
    //return -(2.0-4*pi*pi*x*(x-1.0))*sin(2.0*pi*y);
	return 5.0*pi*pi*sin(pi*x)*sin(2.0*pi*y);
    //return 1.0;
}
void transpose(Matrix bt, Matrix b){
	int i,j,k;
	int len = bt.displ[bt.comm_size-1]+bt.count[bt.comm_size-1];
	double* temp = (double*)malloc(len*sizeof(double));
	double* temp2 = (double*)malloc(len*sizeof(double));
	int l = 0;
	int count = 0;
	for (i = 0; i < bt.comm_size; i++){
		for (j = 0; j < bt.sizes[bt.comm_rank]; j++){
			for (k = 0; k < bt.sizes[i]; k++){
				temp[count++] = b.data[j][k+l];
			}
		}
		l += bt.sizes[i];
	}

	MPI_Alltoallv(temp,b.count,b.displ,MPI_DOUBLE,temp2,bt.count,bt.displ,MPI_DOUBLE,MPI_COMM_WORLD);
	count = 0;
	for (i = 0; i < bt.m; i++){
		for (j = 0; j < bt.sizes[bt.comm_rank]; j++){
			bt.data[j][i] = temp2[count++];
		}
	}
	//free(temp);
	//free(temp2);

}

int main(int argc, char **argv ){

    const double pi = 4.0*atan(1.0);
    double *lambda, **z, *globalDispl;
	Matrix u, ut;
    int rank, size, tag;

    MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc < 2) {
        printf("usage: %s <N> [L]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }

    int N  = atoi(argv[1]);
    int M  = N-1;
    double L=1;
    if (argc > 2)
        L = atof(argv[2]);

    int threads = omp_get_max_threads();
    printf("Threads: %i\n",threads);

    int NN = 4*N;
    double h = L/N;
    int m = M/size; //m = per processor

    u.count = ut.count = calloc(size, sizeof(int));
	u.sizes = ut.sizes = calloc(size, sizeof(int));
	u.displ = ut.displ = calloc(size, sizeof(int));
	u.M = ut.M = M;
	u.comm_size = ut.comm_size = size;
	u.comm_rank = ut.comm_rank = rank;
    u.data = createMatrix(m,M);
    ut.data = createMatrix(m,M);

    globalDispl = createVector(size);
	lambda    = createVector(M);
	z       = createMatrix(threads, NN);

    printf("Eigenvalues\n");
	for (int i = 0; i < M; i++){
		lambda[i] = 2.0*(1.0-cos((i+1)*pi*h));
	}
    printf("G\n");
#pragma omp parallel for schedule(static)
	for (int i = 0; i < m; i++){
		for (int j = 0; j < M; j++){
            //use scaleVec??
			u.data[i][j] = h*h*funcEval(i, j, h, globalDispl[rank], pi);
		}
	}
    // Find out how many rows each processor should have.
    printf("Displacement\n");
	for (int i = 0; i < size; i++){
		if (i < M-m*size){
			u.sizes[i] = ut.sizes[i] = m+1;
		} else {
			u.sizes[i] = ut.sizes[i] = m;
		}
		if (i > 0) {
            globalDispl[i] = globalDispl[i-1] + u.sizes[i-1];
		}
	}

    // Find out how many elements each processor should send/recv in the transpose,
	// and each processor's displacement from the start of the send buffer.
	int sum = 0;
	for (int i = 0; i < size; i++){
		u.count[i] = ut.count[i] = u.sizes[rank]*u.sizes[i];
		u.displ[i] = ut.displ[i] = sum;
		sum = sum + u.count[i];
	}
	// Find out how many rows the current processor gets.
	if ((m*size < M) && (rank < M-m*size)){
		m++;
	}


  //Vector grid = createVector(M);
  //for (int i=0;i<M;++i)
   // grid->data[i] = (i+1)*h;

  //evalMesh(u->as_vec, grid, grid, poisson_source);
  //scaleVector(u->as_vec, h*h);

  //Vector z = createVector(NN);

  //double time = WallTime();

    double startTime = MPI_Wtime();
    printf("Step 1\n");
// Step 1) G = QtGQ = S⁻¹((S(G))t)
#pragma omp parallel for schedule(static)
    for (int i=0; i < m; i++){
        printf("m: %i \t i: %i\n",m,i);
        fst(u.data[i], &N, z[omp_get_thread_num()], &NN);
    }
printf("trans...");
    transpose(ut, u);
printf("posed!\n");
#pragma omp parallel for schedule(static)
    for (int i=0; i < m; i++){
        fstinv(ut.data[i], &N, z[omp_get_thread_num()], &NN);
    }
    printf("Step 2\n");
// Step 2) u_ij = g_ij / (l_i + l_j)
    for (int i=0; i < m; i++){
        for (int j=0; j < M; j++){
            ut.data[i][j] /= lambda[(int)globalDispl[rank]+i]+lambda[j];
        }
    }
    printf("Step 3\n");
// Step 3) U = QUQ^T = S⁻¹((S(Ut))t)
    for (int i=0; i < m; i++)
        fst(ut.data[i], &N, z[omp_get_thread_num()], &NN);

    transpose(u, ut);

    for (int i=0; i < m; i++)
        fstinv(u.data[i], &N, z[omp_get_thread_num()], &NN);

    //evalMesh2(u->as_vec, grid, grid, exact_solution, -1.0);
    //double max = maxNorm(u->as_vec);

    double max = findMax(u);

    if (rank == 0) {
        printf("elapsed: %f\n", MPI_Wtime()-startTime);
        printf("max: %f\n", max);
    }

    MPI_Finalize();
    return 0;
}
