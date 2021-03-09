#define MIN(a,b)  (((a<b)?a:b))


static inline unsigned long long cycles() 
{
  unsigned long long u;
  asm volatile ("rdtscp;shlq $32,%%rdx;orq %%rdx,%%rax;movq %%rax,%0":"=q"(u)::"%rax", "%rdx", "rcx");
  return u;
}

//

double myseconds()
{
    struct timeval tp;
    struct timezone tzp;
    int i;

    i = gettimeofday(&tp,&tzp);
    return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

//

void init(double*  p, int n, int m)
{
#pragma omp parallel for
        for (int j = 0;  j < n; ++j)
        {
                for (int i = 0; i < n; i++)
                {

                        if ( ((i == 0) || (j == 0)) || (i == n-1) || (j == m-1) )
                                p[j*n + i] = 1.;
                        else
                                p[j*n + i] = 0.;
                }
        }
}


double maxNorm(double* v1, double* v2, int size)
{

        double mymax = 0.;
#pragma omp parallel for reduction(max: mymax)
        for (int ii = 0; ii < size; ++ii)
        {
                if (fabs(*v1 - *v2) > mymax)
                {
                        mymax = fabs(*v1 - *v2);
                }
                ++v1; ++v2;
        }
        return mymax;
}


double l2Norm(double* v1, double* v2, int size)
{

        double myl2 = 0.;
#pragma omp parallel for reduction(+: myl2)
        for (int ii = 0; ii < size; ++ii)
        {
                myl2 += fabs(v1[ii]-v2[ii])*(v1[ii] - v2[ii]);
        }
        return sqrt(myl2)/size;
}




void print(double* p, int m, int n)
{
        for (int i=0; i < MIN(n, 15); ++i)
        {
                for (int j=MIN(m, 15); j > 0; --j)
                {
                        printf("%e ", *p);
                        ++p;
                }
                p += m - MIN(m, 15);
                printf("\n");
        }
}
