/*
  This exercise is taken from the class Parallel Programming Workshop (MPI,
  OpenMP and Advanced Topics) at HLRS given by Rolf Rabenseifner
 */

#include <chrono>
#include <cstdio>
#include <cmath>
#include <mpi.h>

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
  int psize, prank;

  MPI_Init(NULL, NULL);

  MPI_Comm_size(MPI_COMM_WORLD, &psize);
  MPI_Comm_rank(MPI_COMM_WORLD, &prank);

  auto mpi_t1 = MPI_Wtime();
  auto t1 = clk::now();

  /* calculate pi = integral [0..1] 4 / (1 + x**2) dx */
  dx = 1. / n;
  sum = 0.0;
  for (i = 1; i <= n; i++) {
    x = (1. * i - 0.5) * dx;
    sum = sum + f(x);
  }
  pi = dx * sum;



  auto mpi_elapsed = MPI_Wtime() - mpi_t1;
  second elapsed = clk::now() - t1;

  if(prank == 0) {
    std::printf("computed pi                 = %.16g\n", pi);
    std::printf("wall clock time (mpi_wtime) = %.4gs with %d process\n", mpi_elapsed, psize);
    std::printf("wall clock time (chrono)    = %.4gs\n", elapsed.count());

    for(int d = 1; d <= 15; ++d) {
      std::printf("%d", digit(pi, d));
    }
  }

  MPI_Finalize();

  return 0;
}
