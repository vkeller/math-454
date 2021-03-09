#include "io_binary_mpi_io.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#define epsilon 0.005

static inline void compute_row(int i, int n, float ** uo, float ** u,
			       float ** f, float * l2, float h) {
  int j;
  for (j = 1; j < n; j++) {
    // computation of the new step
    u[i][j] = 0.25 * (uo[i - 1][j] + uo[i + 1][j] + uo[i][j - 1] +
		      uo[i][j + 1] - f[i][j] * h * h);

    // L2 norm
    *l2 += (uo[i][j] - u[i][j]) * (uo[i][j] - u[i][j]);
  }
}


static inline void swap_grids(int n, int n_loc, float ** uo, float ** u) {
  int i, j;

  for (i = 0; i < n_loc; i++) {
    for (j = 0; j < n + 1; j++) {
      uo[i][j] = u[i][j];
    }
  }
}

enum _requests_t {
  _north_send = 0,
  _sourth_send = 1,
  _north_recv = 2,
  _sourth_recv = 3,
  _max_requests_t = 4,
};

int main(int argc, char * argv[]) {
  int i, j, k = 0;
  float l2;
  float **u, **uo, **f;
  float *u_data, *uo_data, *f_data;
  float h;

  int prank, psize;
  int n, n_loc, offset;

  MPI_Status status[4];
  MPI_Request requests[4];

  int N = 1024;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &prank);
  MPI_Comm_size(MPI_COMM_WORLD, &psize);

  if(argc == 2) {
    N = atoi(argv[1]);
  }

  h = 1. / (N + 1);

  // n is N + 1 minus the first and last line that are forced to 0 in the
  // computation. This 2 lines will be replaced by ghosts lines
  n = N - 1;

  // Divide n in psize and distribute the excess to the (n % psize) proc
  n_loc = (n / psize) + (prank < n % psize ? 1 : 0);

  {
    // Computing the offset of where in the global array the local array is
    // located. (this is needed to initialize f) It could be computed locally
    // without communication but this way you see a Allgather

    // count the first line in the offset to compute correctly f
    int n_loc_tmp = n_loc + (prank == 0 ? 1 : 0);

    int buf[psize];
    MPI_Allgather(&n_loc_tmp, 1, MPI_INT, buf, 1, MPI_INT, MPI_COMM_WORLD);

    offset = 0;
    for(i = 0; i < prank; ++i) {
      offset += buf[i];
    }
  }

  // North and south ghost lines
  n_loc += 2;

  // Allocation of the data storage
  u_data  = (float *) malloc((N + 1) * n_loc * sizeof(float));
  uo_data = (float *) malloc((N + 1) * n_loc * sizeof(float));
  f_data  = (float *) malloc((N + 1) * n_loc * sizeof(float));

  // Allocation of arrays of pointers for rows
  u  = (float **) malloc(n_loc * sizeof(float *));
  uo = (float **) malloc(n_loc * sizeof(float *));
  f  = (float **) malloc(n_loc * sizeof(float *));

  // set the row pointers in the memory
  for (i = 0; i < n_loc; i++) {
    u[i]  = u_data  + i * (N + 1);
    uo[i] = uo_data + i * (N + 1);
    f[i]  = f_data  + i * (N + 1);
  }

  int i_start = prank == 0 ? 0 : 1;
  int i_end   = prank == (psize - 1) ? n_loc : n_loc -1;

  // initialization of u0 and f
  for (i = i_start; i < i_end; i++) {
    for (j = 0; j < N + 1; j++) {
      float x = (i - i_start + offset) * h;
      float y = (    j) * h;
      u[i][j] = 0;
      uo[i][j] = 0;
      f[i][j] = -2. * 100. * M_PI * M_PI * sin(10. * M_PI * x) *
	sin(10. * M_PI * y);
    }
  }

  double start = MPI_Wtime();

  int north_proc = prank == 0 ? MPI_PROC_NULL : prank - 1;
  int south_proc = prank == (psize - 1) ? MPI_PROC_NULL : prank + 1;

  // initializing the persistent communications
  /// synchronize north
  MPI_Send_init(uo[1], N + 1, MPI_FLOAT, north_proc, 0, MPI_COMM_WORLD,
		requests + _north_send);
  MPI_Recv_init(uo[0], N + 1, MPI_FLOAT, north_proc, 0, MPI_COMM_WORLD,
		requests + _north_recv);

  /// synchronize south
  MPI_Send_init(uo[n_loc - 2], N + 1, MPI_FLOAT, south_proc, 0, MPI_COMM_WORLD,
		requests + _sourth_send);
  MPI_Recv_init(uo[n_loc - 1], N + 1, MPI_FLOAT, south_proc, 0, MPI_COMM_WORLD,
		requests + _sourth_recv);

  do {
    l2 = 0.;

    MPI_Startall(_max_requests_t, requests);

    /// do not treat the layers adjacent to ghost regions
    for (i = 2; i < n_loc - 2; i++) {
      compute_row(i, N, uo, u, f, &l2, h);
    }

    /// wait to receive the ghosts before using them for them computation
    MPI_Waitall(2, requests + _north_recv, status);

    /// treat the layers adjacent to ghost regions
    compute_row(1, N, uo, u, f, &l2, h);
    compute_row(n_loc - 2, N, uo, u, f, &l2, h);

    // reduce the 12 norm to every proc
    MPI_Allreduce(MPI_IN_PLACE, &l2, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

    /// wait to send everything before changing the buffers
    MPI_Waitall(2, requests + _north_send, status);

    // copy new grid in old grid
    swap_grids(N, n_loc, uo, u);

//    write_to_file(N + 1,
//		  i_end - i_start, // to remove the ghosts
//		  u_data + i_start * (N + 1), // to start at the first not ghost
//		  // line moves the pointer of
//		  // i_start lines
//		  k, -1., 1.);

    ++k;
  } while (l2 > epsilon);

  double ttime = MPI_Wtime() - start;

  if(prank == 0) {
    printf("Execution time = %f [s] (%f [ms/step]) - nb steps %d [l2 = %g] - nb procs %d\n",
	   ttime, (ttime / k) * 1e3, k, l2, psize);
  }

  for(i = 0; i < 4; ++i) {
    MPI_Request_free(requests + i);
  }

  // de-allocate the data
  free(u_data);
  free(uo_data);
  free(f_data);

  // de-allocate the rows pointers
  free(u);
  free(uo);
  free(f);

  MPI_Finalize();

  return 0;
}
