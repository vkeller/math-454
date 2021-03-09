#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <mpi.h>

#define N 255

#define NORTH 0
#define SOUTH 1

#define eps 0.005

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

int main(int argc, char *argv[]) {
  int     i, j, k;
  float **u;
  float **uo;
  float **f;
  float   h  = 0.;
  float   l2 = 0.;

  MPI_Status status;
  double start, end;

  int n, n_loc, offset;
  int prank, psize;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &prank);
  MPI_Comm_size(MPI_COMM_WORLD, &psize);

  n = N + 1;
  h = 1. / n;

  // Divide n in psize and distribute the excess to the (n % psize) proc
  n_loc = (n / psize) + (prank < n % psize ? 1 : 0);

  { // Computing the offset of where in the global array the local array is
    // located. (this is needed to initialize f) It could be computed locally
    // without communication but this way you see a Allgather
    int i;
    int buf[psize];
    MPI_Allgather(&n_loc, 1, MPI_INT, buf, 1, MPI_INT, MPI_COMM_WORLD);

    offset = 0;
    for(i = 0; i < prank; ++i) {
      offset += buf[i];
    }
  }

  // add 2 for north and south ghost
  n_loc += 2;


  // Allocation of arrays of pointers to the number of rows per proc
  u  = (float**) malloc(n_loc * sizeof(float*));
  uo = (float**) malloc(n_loc * sizeof(float*));
  f  = (float**) malloc(n_loc * sizeof(float*));

  // allocation of the rows to the size of full lines
  for(i=0; i < n_loc; i++) {
    u [i] = (float*) malloc(n*sizeof(float));
    uo[i] = (float*) malloc(n*sizeof(float));
    f [i] = (float*) malloc(n*sizeof(float));
  }


  // initialization of u0 and f
  for(i = 1; i < n_loc - 1; i++) {
    for(j = 0; j < n; j++) {
      u [i][j] = 0;
      uo[i][j] = 0;
      f [i][j] = -2.*100. * M_PI * M_PI * sin(10.*M_PI*((i-1) + offset)*h) * sin(10.*M_PI*j*h);
    }
  }

  k=0;

  start = MPI_Wtime();
  do {
    int i_start = 1;
    int i_end   = n_loc - 1;

    l2 = 0.;

    // First synchronize uo
    if(prank > 0) { // send recv with top proc
      MPI_Sendrecv(uo[1], n, MPI_FLOAT, prank - 1, NORTH,
                   uo[0], n, MPI_FLOAT, prank - 1, SOUTH, MPI_COMM_WORLD, &status);
    } else {
      ++i_start;
    }

    if(prank < psize - 1) { // send recv with bottom proc
      MPI_Sendrecv(uo[n_loc - 2], n, MPI_FLOAT, prank + 1, SOUTH,
                   uo[n_loc - 1], n, MPI_FLOAT, prank + 1, NORTH, MPI_COMM_WORLD, &status);
    } else {
      --i_end;
    }

    for(i = i_start; i < i_end; i++) {
      for(j = 1; j < n - 1; j++) {
        // computation of the new step
        u[i][j] = 0.25 * ( uo[i-1][j] + uo[i+1][j] + uo[i][j-1] + uo[i][j+1] - f[i][j]*h*h);

        // L2 norm
        l2 += (uo[i][j] - u[i][j])*(uo[i][j] - u[i][j]);
      }
    }

    // copy new grid in old grid
    for(i = 1; i < n_loc - 1; i++) {
      for(j = 0; j < n; j++){
        uo[i][j] = u[i][j];
      }
    }


    // reduce the 12 norm to every proc
    MPI_Allreduce(MPI_IN_PLACE, &l2, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

    // outputs
    if(prank == 0) {
      printf("l2=%.5f (k=%d)\n", sqrt(l2), k);
    }

    k++;
  } while(l2 > eps);
  end = MPI_Wtime();

  if(prank == 0) {
    double t = end - start;
    printf("T=%.5f s (%.5f s/step)\n", t, t/k);
  }

  // deallocation of the rows
  for(i=0; i < n_loc ;i++) {
    free(u [i]);
    free(uo[i]);
    free(f [i]);
  }

  // deallocate the pointers
  free(u);
  free(uo);
  free(f);

  MPI_Finalize();

  return 0;
}
