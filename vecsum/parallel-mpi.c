#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <sys/time.h>

double walltime ( void ) {  // Returns current time, with microsecond precision
#ifdef HAVE_MPI
    return MPI_Wtime();
#endif
    static struct timeval t;
    gettimeofday ( &t, NULL );
    return ( t.tv_sec + 1e-6 * t.tv_usec );
}

int isPowerOfTwo (unsigned int x)
{
  return ((x != 0) && !(x & (x - 1)));
}

double time_init;

int main(int argc, char** argv)
{
    double pi = 4.0*atan(1);
    double sum = pi*pi/6;
    int rank = 0,size = 1;
    
    if (argc < 2) {
        printf("Need one parameter, the size of the vector\n");
        return 1;
    }
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int N = atoi(argv[1]);
    int share = N/size;
    double sumn = 0.0;
    double* v = (double*)calloc(share,sizeof(double));

    if (rank == 0){
        printf("MPI   \tThreadcount: %i\n",size);
    }
    if (!isPowerOfTwo(size)){
        if (rank == 0){
        printf("The number of processors must be a power of 2, %i processors supplied.\n",size);
        }
        MPI_Finalize();
        return 1;
    }
    time_init = walltime();

    for (long int i=share; i>0; i--){
        v[i-1] = (double)1.0/((i+rank*share)*(i+rank*share));
        sumn += v[i-1];
    }

    double s2 = sumn;
    MPI_Reduce(&s2, &sumn, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);

    if (rank == 0){
        printf("Error:\t\t%.52f\nTime Elapsed:\t%f\n",sum-sumn,walltime()-time_init);
    }
    MPI_Finalize();

return 0;
}


