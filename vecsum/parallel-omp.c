#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>

double walltime ( void ) {  // Returns current time, with microsecond precision
    static struct timeval t;
    gettimeofday ( &t, NULL );
    return ( t.tv_sec + 1e-6 * t.tv_usec );
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

    double* v = (double*)malloc(n*sizeof(double));
    double sumn = 0;

    time_init = walltime();

#pragma omp parallel for schedule(static) reduction(+:sumn)
    for(int i=1; i<n+1; i++){
        v[i-1] = (double)1.0/(i*i);
        sumn += v[i-1];
    }
    printf("Error:\t\t%.52f\nTime Elapsed:\t%f\n",sum-sumn,walltime()-time_init);

return 0;
}
/*
#pragma omp parallel for schedule(static) reduction(+:sumn)
    for(int i=1; i<(n/P +1); i++){
        v[omp_get_thread_num()*n/P + i-1] = (double)1.0/(i*i);
        sumn += v[omp_get_thread_num()*n/P + i-1];
        printf("Thread: %i\tIndex: %i\n",omp_get_thread_num(),omp_get_thread_num()*n/P + i-1);
    }
*/
