#include <sys/time.h>

#define N 2056
#define eps 0.005

//#define _OUTPUT

double start, end;

double second()
{
        struct timeval tp;
        struct timezone tzp;
        gettimeofday(&tp,&tzp);
        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}


