#include "timing.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

long int fibonacci(int n) {
  long int x, y;
  if (n < 2) {
    return n;
  } else {
#pragma omp task shared(x)
    x = fibonacci(n-1);
#pragma omp task shared(y)
    y = fibonacci(n-2);
#pragma omp taskwait
    return (x+y);
  }
}

int main() {
  int n = 42;
  double t1, t2;
  long int fib = 0;
  
  t1 = second();
#pragma omp parallel
  {
#pragma omp single nowait
    {
      fib = fibonacci(n);
    }
  }
  
  t2 = second();

  printf("fib(%d) = %ld (in %g [s])\n", n, fib, (t2-t1));
  
  return 0;
}
