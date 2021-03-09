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
  double dx, x, pi;
  int psize, prank;

  struct Sum {
    double sum;
    int rank;
  };

  MPI_Init(NULL, NULL);

  MPI_Comm_size(MPI_COMM_WORLD, &psize);
  MPI_Comm_rank(MPI_COMM_WORLD, &prank);

  Sum sum{0., 0};
  Sum send, recv;

  int blk_length[2] = {1, 1};

  MPI_Aint zero_address, first_address, second_address;
  MPI_Get_address(&send, &zero_address);
  MPI_Get_address(&send.sum, &first_address);
  MPI_Get_address(&send.rank, &second_address);

  MPI_Aint displs[2];
  displs[0] = MPI_Aint_diff(first_address, zero_address);;
  displs[1] = MPI_Aint_diff(second_address, first_address);

  MPI_Datatype types[2] = {MPI_DOUBLE, MPI_INT};
  MPI_Datatype sum_t;
  MPI_Type_create_struct(2, blk_length, displs, types, &sum_t);
  MPI_Type_commit(&sum_t);

  auto mpi_t1 = MPI_Wtime();
  auto t1 = clk::now();

  int nlocal = n / psize;
  int istart = 1 + nlocal * prank;
  int iend = nlocal * (prank + 1);

  /* calculate pi = integral [0..1] 4 / (1 + x**2) dx */
  dx = 1. / n;

  for (i = istart; i <= iend; i++) {
    x = (1. * i - 0.5) * dx;
    sum.sum = sum.sum + f(x);
  }

  int next = (prank + 1) % psize;
  int prev = (prank - 1 + psize) % psize;

  MPI_Request request;
  send = sum;
  send.rank = prank;
  for(int s = 1; s < psize; ++s) {
    MPI_Isend(&send, 1, sum_t, next, 13, MPI_COMM_WORLD, &request);
    MPI_Recv(&recv, 1, sum_t, prev, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    sum.sum += recv.sum;
    sum.rank += recv.rank;

    MPI_Wait(&request, MPI_STATUS_IGNORE);

    send = recv;
  }

  pi = dx * sum.sum;

  auto mpi_elapsed = MPI_Wtime() - mpi_t1;
  second elapsed = clk::now() - t1;

  if(prank == 0) {
    std::printf("computed pi                 = %.16g / %d\n", pi, sum.rank);
    std::printf("wall clock time (mpi_wtime) = %.4gs with %d process\n", mpi_elapsed, psize);
    std::printf("wall clock time (chrono)    = %.4gs\n", elapsed.count());

    for(int d = 1; d <= 15; ++d) {
      std::printf("%d", digit(pi, d));
    }
  }

  MPI_Type_free(&sum_t);

  MPI_Finalize();

  return 0;
}
