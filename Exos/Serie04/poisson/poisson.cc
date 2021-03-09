/* -------------------------------------------------------------------------- */
#include "simulation.hh"
/* -------------------------------------------------------------------------- */
#include <chrono>
#include <iostream>
#include <tuple>
/* -------------------------------------------------------------------------- */

#define N 256
#define EPSILON 0.005

int main() {

  Simulation simu(N, N);

  simu.set_initial_conditions();

  simu.set_epsilon(EPSILON);

  float l2;
  int k;
  std::tie(l2, k) = simu.compute();

  std::cout <<"nb steps " << k << " [l2 = " << l2 << "]"
            << std::endl;

  return 0;
}
