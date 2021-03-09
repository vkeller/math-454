/*
  This exercise is taken from the class Parallel Programming Workshop (MPI,
  OpenMP and Advanced Topics) at HLRS given by Rolf Rabenseifner
 */

#include <chrono>
#include <cstdio>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

using clk = std::chrono::high_resolution_clock;
using second = std::chrono::duration<double>;
using time_point = std::chrono::time_point<clk>;

inline int digit(double x, int n) {
  return std::trunc(x * std::pow(10., n)) - std::trunc(x * std::pow(10., n - 1)) *10.;
}

inline double f(double a) { return (4. / (1. + a * a)); }

const int n = 10000000;

int main(int /* argc */ , char ** /* argv */) {
  int i;
  double dx, x, sum, pi;
  int nthreads;

#ifdef _OPENMP
  nthreads = omp_get_max_threads();
  auto omp_t1 = omp_get_wtime();
#endif
  auto t1 = clk::now();

  /* calculate pi = integral [0..1] 4 / (1 + x**2) dx */
  dx = 1. / n;
  sum = 0.0;
#pragma omp parallel for
  for (i = 0; i < n; i++) {
    x = 1. * i * dx;
    sum = sum + f(x);
  }
  pi = dx * sum;



#ifdef _OPENMP
  auto omp_elapsed = omp_get_wtime() - omp_t1;
#endif
  second elapsed = clk::now() - t1;


  std::printf("computed pi                     = %.16g\n", pi);
#ifdef _OPENMP
  std::printf("wall clock time (omp_get_wtime) = %.4gs in %d threads\n", omp_elapsed, nthreads);
#endif
  std::printf("wall clock time (chrono)        = %.4gs\n", elapsed.count());

  for(int d = 1; d <= 15; ++d) {
    std::printf("%d", digit(pi, d));
  }

  return 0;
}
