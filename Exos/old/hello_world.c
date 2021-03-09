#include <stdio.h>
#include <omp.h>
int main(int argc, char *argv[]) {
   int myrank=0;
   int mysize=1;
#if defined (_OPENMP)
#pragma omp parallel default(shared) private(myrank, mysize)
{
   mysize = omp_get_num_threads();
   myrank = omp_get_thread_num();
#endif
   printf("Hello from thread %d out of %d\n", myrank, mysize);
#if defined (_OPENMP)
}
#endif
   return 0;
}
