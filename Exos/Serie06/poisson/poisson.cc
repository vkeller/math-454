/* -------------------------------------------------------------------------- */
#include "simulation.hh"
/* -------------------------------------------------------------------------- */
#include <chrono>
#include <iostream>
#include <sstream>
#include <tuple>
#include <chrono>
/* -------------------------------------------------------------------------- */
#include <omp.h>
#include <mpi.h>
/* -------------------------------------------------------------------------- */

#define EPSILON 0.005

typedef std::chrono::high_resolution_clock clk;
typedef std::chrono::duration<double> second;

static void usage(const std::string & prog_name) {
  std::cerr << prog_name << " <grid_size>" << std::endl;
  exit(0);
}


int main(int argc, char * argv[]) {
  MPI_Init(&argc, &argv);
  int prank, psize;

  MPI_Comm_rank(MPI_COMM_WORLD, &prank);
  MPI_Comm_size(MPI_COMM_WORLD, &psize);

  if (argc != 2) usage(argv[0]);

  std::stringstream args(argv[1]);
  int N;
  args >> N;

  if(args.fail()) usage(argv[0]);

  Simulation simu(N, N, MPI_COMM_WORLD);

  simu.set_initial_conditions();

  simu.set_epsilon(EPSILON);

  float l2;
  int k;

  auto start = clk::now();
  std::tie(l2, k) = simu.compute();
  auto end = clk::now();

  second time = end - start;

  if(prank == 0)
    std::cout << psize << " " << N << " "
              << k << " " << std::scientific << l2 << " "
              << time.count() << std::endl;

  MPI_Finalize();

  return 0;
}
