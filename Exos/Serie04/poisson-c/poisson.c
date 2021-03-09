/* -------------------------------------------------------------------------- */
#include "timing.h"
#include "simulation.h"
/* -------------------------------------------------------------------------- */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/* -------------------------------------------------------------------------- */

#define epsilon 0.005

int main(int argc, char * argv[]) {
  int k = 0;
  float l2;
  float **u, **uo, **f;
  float h;

  int n = 256;

  if (argc == 2) {
    n = atoi(argv[1]);
  }

  h = 1. / n;

  allocate_grids(n, &uo, &u, &f);

  // initialization of u0 and f
  initialize_grids(n, uo, u, f, h);
  
  double start = second();

  l2 = simulate(n, uo, u, f, h, epsilon, &k);
  
  double ttime = second() - start;
  printf("%d %d %f %f\n", n, k, l2, ttime);

  deallocate_grids(&uo, &u, &f);

  return 0;
}
