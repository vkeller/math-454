#include "timing.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

long int fibonacci_seq(int n) {
  long int x, y;
  if (n < 2) {
    return n;
  } else {
    x = fibonacci_seq(n-1);
    y = fibonacci_seq(n-2);
    return (x+y);
  }
}

long int fibonacci(int n, int level, int cutoff) {
  long int x, y;
  if (n < 2) {
    return n;
  } else if (level < cutoff) {
#pragma omp task shared(x)
    x = fibonacci(n-1, level+1, cutoff);
#pragma omp task shared(y)
    y = fibonacci(n-2, level+1, cutoff);
#pragma omp taskwait
    return (x+y);
  } else {
    x = fibonacci_seq(n - 1);
    y = fibonacci_seq(n - 2);
    return (x+y);
  }
}


int main() {
  int n = 42;
  int cutoff = 10;
  double t1, t2;
  long int fib = 0;
  
  t1 = second();
#pragma omp parallel
  {
#pragma omp single nowait
    {
      fib = fibonacci(n, 0, cutoff);
    }
  }
  
  t2 = second();

  printf("Fib(%d) = %ld (in %g [s])\n", n, fib, (t2-t1));
  
  return 0;
}
