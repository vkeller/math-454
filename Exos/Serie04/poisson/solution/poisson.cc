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
/* -------------------------------------------------------------------------- */

#define EPSILON 0.005

typedef std::chrono::high_resolution_clock clk;
typedef std::chrono::duration<double> second;

static void usage(const std::string & prog_name) {
  std::cerr << prog_name << " <grid_size>" << std::endl;
  exit(0);
}

int main(int argc, char * argv[]) {
  if (argc != 2) usage(argv[0]);

  std::stringstream args(argv[1]);
  int N;
  args >> N;

  if(args.fail()) usage(argv[0]);

  Simulation simu(N, N);

  simu.set_initial_conditions();

  simu.set_epsilon(EPSILON);

  float l2;
  int k;

  auto start = clk::now();
  std::tie(l2, k) = simu.compute();
  auto end = clk::now();

  second time = end - start;

  std::cout << omp_get_max_threads() << " " << N << " "
            << k << " " << std::scientific << l2 << " "
            << time.count() << std::endl;

  return 0;
}
