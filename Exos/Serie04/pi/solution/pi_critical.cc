#include <chrono>
#include <cmath>
#include <cstdio>

#if defined(_OPENMP)
#include <omp.h>
#endif

using clk = std::chrono::high_resolution_clock;
using second = std::chrono::duration<double>;
using time_point = std::chrono::time_point<clk>;

inline int digit(double x, int n) {
  return std::trunc(x * std::pow(10., n)) -
         std::trunc(x * std::pow(10., n - 1)) * 10.;
}

inline double f(double a) { return (4. / (1. + a * a)); }

const int n = 10000000;

int main(int /* argc */, char ** /* argv */) {
  int i;
  double dx, x, sum, pi;

#if defined(_OPENMP)
  int num_threads = omp_get_max_threads();
#endif

#if defined(_OPENMP)
  auto omp_t1 = omp_get_wtime();
#endif
  auto t1 = clk::now();

  sum = 0.;
#pragma omp parallel shared(sum) private(x)
  {
    /* calculate pi = integral [0..1] 4 / (1 + x**2) dx */
    dx = 1. / n;
#pragma omp for
    for (i = 1; i <= n; i++) {
      x = (1. * i - 0.5) * dx;
#pragma omp critical
      sum = sum + f(x);
    }
  }

  pi = dx * sum;

#if defined(_OPENMP)
  auto omp_elapsed = omp_get_wtime() - omp_t1;
#endif
  second elapsed = clk::now() - t1;

  std::printf("computed pi                     = %.16g\n", pi);
#if defined(_OPENMP)
  std::printf("wall clock time (omp_get_wtime) = %.4gs on %d threads\n",
              omp_elapsed, num_threads);
#endif
  std::printf("wall clock time (chrono)        = %.4gs\n", elapsed.count());

  for (int d = 1; d <= 15; ++d) {
    std::printf("%d", digit(pi, d));
  }

  return 0;
}
