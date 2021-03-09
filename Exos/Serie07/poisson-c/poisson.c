/* -------------------------------------------------------------------------- */
#include "timing.h"
#include "simulation.h"
/* -------------------------------------------------------------------------- */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/* -------------------------------------------------------------------------- */
#include <mpi.h>
/* -------------------------------------------------------------------------- */

#define epsilon 0.005

static void usage(char * prog_name) {
  printf("%s <grid_size>\n", prog_name);
  exit(0);
}

int main(int argc, char * argv[]) {
  int k = 0;
  float l2;
  float **u, **uo, **f;
  

  MPI_Init(&argc, &argv);

  int prank, psize;
  MPI_Comm_rank(MPI_COMM_WORLD, &prank);
  MPI_Comm_size(MPI_COMM_WORLD, &psize);

  if (argc != 2) usage(argv[0]);

  int N = atoi(argv[1]);
  float h = 1. / N;

  int m = N / psize + (prank < N % psize ? 1 : 0);
  int n = N;

  // adding the ghosts lines if needed
  if (psize > 1)
    m += (prank == 0 || prank == psize - 1) ? 1 : 2;

  int offset_m =
      (N / psize) * prank + (prank < N % psize ? prank : N % psize);
  int offset_n = 0;

  
  allocate_grids(m, n, &uo, &u, &f);

  // initialization of u0 and f
  initialize_grids(m, n, offset_m, offset_n, uo, u, f, h);
  
  double start = second();

  l2 = simulate(m, n, uo, u, f, h, epsilon, &k);
  
  double ttime = second() - start;
  printf("%d %d %f %f\n", n, k, l2, ttime);

  deallocate_grids(&uo, &u, &f);

  MPI_Finalize();
  
  return 0;
}
