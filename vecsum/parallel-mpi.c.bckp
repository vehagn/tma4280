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

    if (argc < 2) {
        printf("Need one parameter, the size of the vector\n");
        return 1;
    }
    int n = atoi(argv[1]);

    int rank = 0, size = 1;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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

    double* v = (double*)malloc(n*sizeof(double));
    double sumn = 0;
    /*double* sumn = (double*)malloc(size*sizeof(double));
    for (int i = 0; i < size; i++){
        sumn[i] = 0;
    }*/

    time_init = walltime();

    for (int i=n-(size-(rank+1))*n/size; i>rank*n/size; i--){      
	v[i] = (double)1.0/(i*i);
        sumn += v[i];
       // printf("Rank: %i\tIndex: %i\tSumn: %f\n",rank,i,sumn);
    }
    double s2 = sumn; //What happens here?
    MPI_Allreduce(&s2, &sumn, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0){
        /*for (int i = 1; i < size; i++){
            printf("sumn: %f\n",sumn[0]);
            sumn[0] += sumn[i];
        }*/
        printf("Error:\t\t%.52f\nTime Elapsed:\t%f\n",sum-sumn,walltime()-time_init);
    }
    MPI_Finalize();

return 0;
}


