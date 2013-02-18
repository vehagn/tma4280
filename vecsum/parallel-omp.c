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

    printf("OpenMP\tThreadcount: %i\n",omp_get_max_threads());

    double* v = (double*)malloc(n*sizeof(double));
    double sumn = 0;

    time_init = walltime();

#pragma omp parallel for schedule(static) reduction(+:sumn)
    for(long long int i=n; i>0; i--){
        v[i] = 1.0/((double)i*i);
        sumn += v[i];
    }
    printf("Error:\t\t%.52f\nTime Elapsed:\t%f\n",sum-sumn,walltime()-time_init);

return 0;
}
