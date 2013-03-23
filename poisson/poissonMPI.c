#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include "mpi.h"
#include "omp.h"

typedef double Real;

typedef struct {
	int perProc;
	int m;
	Real** data;
	int comm_size;
	int comm_rank;
	int* sizes;
	int* displ;
	int* count;
} Matrix;

/* function prototypes */
Real *createRealArray (int n);
Real **createReal2DArray (int m, int n);
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);
void printToFile(Real **b, int n, int m, int rank);
void transpose(Matrix bt, Matrix b);
void printArray(Real** array,int n,int m);
Real funcEval(int i, int j, Real h, int displ, Real pi);

int main(int argc, char **argv )
{
	if( argc < 2 ) {
		printf("need a problem size\n");
		return 1;
	}
	
	double start = MPI_Wtime();

	Real *diag, **z, *globalDispl;
	Real pi, h, umax;
	Matrix b, bt;
	int i, j, n, m, nn;
	int rank, nproc, tag;
	n  = atoi(argv[1]);
	m  = n-1;
	nn = 4*n;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// Initialize the number of threads available to the program.
	int threads;
        #pragma omp parallel
        {
                threads = omp_get_num_threads();
        }
	b.count = bt.count = calloc(nproc, sizeof(int));
	b.sizes = bt.sizes = calloc(nproc, sizeof(int));
	b.displ = bt.displ = calloc(nproc, sizeof(int));
	b.m = bt.m = m;
	b.comm_size = bt.comm_size = nproc;
	b.comm_rank = bt.comm_rank = rank;
	
	globalDispl = createRealArray(nproc);
	diag    = createRealArray(m);
	z       = createReal2DArray(threads, nn);
	h       = 1./(Real)n;
	pi      = 4.*atan(1.);
	
	// Find out how many rows each processor should have.
	int perProc = m/nproc;
	for (i = 0; i < nproc; i++){
		if (i < m-perProc*nproc){
			b.sizes[i] = bt.sizes[i] = perProc+1;
		} else {
			b.sizes[i] = bt.sizes[i] = perProc;
		}
		if (i > 0) globalDispl[i] = globalDispl[i-1] + b.sizes[i-1];
	}
	// Find out how many elements each processor should send/recv in the transpose,
	// and each processor's displacement from the start of the send buffer.
	int sum = 0;	
	for (i = 0; i < nproc; i++){
		b.count[i] = bt.count[i] = b.sizes[rank]*b.sizes[i];
		b.displ[i] = bt.displ[i] = sum;
		sum = sum + b.count[i];
	}
	// Find out how many rows the current processor gets.	
	if ((perProc*nproc < m) && (rank < m-perProc*nproc)){
		perProc++;
	}

	b.data  = createReal2DArray(perProc,m);
	bt.data = createReal2DArray(perProc,m);

	for (i = 0; i < m; i++){
		diag[i] = 2.*(1.-cos((i+1)*pi*h));
	}
	
	Real h2 = h*h;
	#pragma omp parallel for schedule(static)
	for (i = 0; i < perProc; i++){
		for (j = 0; j < m; j++){
			b.data[i][j] = h2*funcEval(i, j, h, globalDispl[rank], pi);
		}
	}
if (rank == 0){	
printf("Init\n");
printf("data: %f\n",b.data[0][0]);
}
	#pragma omp parallel for schedule(static)
	for (i = 0; i < perProc; i++){
		fst_(b.data[i], &n, z[omp_get_thread_num()], &nn);
	}

	transpose(bt,b);
	
	#pragma omp parallel for schedule(static)
	for (i = 0; i < perProc; i++) {
		fstinv_(bt.data[i], &n, z[omp_get_thread_num()], &nn);
	}
if (rank == 0){
printf("Step 1\n");
printf("data: %f\n",bt.data[0][0]);
}
	
	#pragma omp parallel for schedule(static)
	for (i = 0; i < perProc; i++) {
		for (j = 0; j < m; j++) {
			bt.data[i][j] = bt.data[i][j]/(diag[(int)globalDispl[rank]+i]+diag[j]);
		}
	}
if (rank == 0){
printf("Step 2\n");
printf("data: %f\n",bt.data[0][0]);
}
	#pragma omp parallel for schedule(static)
	for (i=0; i < perProc; i++) {
		fst_(bt.data[i], &n, z[omp_get_thread_num()], &nn);
	}
	
	transpose (b,bt);

	#pragma omp parallel for schedule(static)
	for (i = 0; i < perProc; i++) {
		fstinv_(b.data[i], &n, z[omp_get_thread_num()], &nn);
	}

	umax = 0.0;
	for (i = 0; i < perProc; i++) {
	 	for (j = 0; j < m; j++) {
	 		if (b.data[i][j]-funcEval(i, j, h, globalDispl[rank], pi)/(5.*pi*pi) > umax)
				 umax = b.data[i][j]-funcEval(i, j, h, globalDispl[rank], pi)/(5.*pi*pi);
	 	}
	}
	MPI_Allreduce(&umax, &umax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);	
	//printf (" Max pointwise error = %.14f \n", umax);

	printToFile(b.data, perProc, m, rank);	
	//printf("Rank: %d of %d, perProc: %d, time: %.15f\n", rank, nproc, perProc, MPI_Wtime()-start);


	MPI_Finalize();
	return 0;
}

void printToFile(Real **b, int n, int m, int rank){
	FILE *file;
	char filename[20];
	sprintf(filename, "out%d.dat", rank);
	file = fopen(filename, "w+");
	int i,j;
	for (i = 0; i < n; i++){
		for (j = 0; j < m; j++){
			fprintf(file, "%.10f ", b[i][j]);
		}
		fprintf(file, "\n");
	}
}

Real *createRealArray (int n)
{
	Real *a;
	int i;
	a = (Real *)malloc(n*sizeof(Real));
	for (i=0; i < n; i++) {
		a[i] = 0.0;
	}
	return (a);
}

Real **createReal2DArray (int n1, int n2)
{
	int i, n;
	Real **a;
	a    = (Real **)malloc(n1   *sizeof(Real *));
	a[0] = (Real  *)malloc(n1*n2*sizeof(Real));
	for (i=1; i < n1; i++) {
		a[i] = a[i-1] + n2;
	}
	n = n1*n2;
	memset(a[0],0,n*sizeof(Real));
	return (a);
}

void transpose(Matrix bt, Matrix b)
{
	int i,j,k;
	int len = bt.displ[bt.comm_size-1]+bt.count[bt.comm_size-1];
	Real* temp = (Real*)malloc(len*sizeof(Real));
	Real* temp2 = (Real*)malloc(len*sizeof(Real));
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
}

void printArray(Real** array, int n, int m)
{
	int i,j;
	for (i=0; i<n ; i++) {
		for (j = 0; j < m; j++) {
			printf("%f ",array[i][j]);
		}
		printf("\n");
	}
}

Real funcEval(int i, int j, Real h, int displ, Real pi)
{
	Real y = (Real)(displ+i+1)*h;
	Real x = (Real)(j+1)*h;
//	if (x > 0.6 && y > 0.6) return 10.;
//	return exp(1.-x-y)-1.;		
	return 5.*pi*pi*sin(pi*x)*sin(2.*pi*y);
//	return 1.0;
}
