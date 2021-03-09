/* -------------------------------------------------------------------------- */
#include "io_binary.h"
/* -------------------------------------------------------------------------- */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/* -------------------------------------------------------------------------- */
#include <mpi.h>
/* -------------------------------------------------------------------------- */

void allocate_grids(int m, int n, float ***u, float ***uo, float ***f) {
  int i;

  // Allocation of the data storage
  float *u_data = (float *)malloc(m * n * sizeof(float));
  float *uo_data = (float *)malloc(m * n * sizeof(float));
  float *f_data = (float *)malloc(m * n * sizeof(float));

  // Allocation of arrays of pointers for rows
  *u = (float **)malloc(m * sizeof(float *));
  *uo = (float **)malloc(m * sizeof(float *));
  *f = (float **)malloc(m * sizeof(float *));

  // set the row pointers in the memory
  for (i = 0; i < m; i++) {
    (*u)[i] = u_data + i * n;
    (*uo)[i] = uo_data + i * n;
    (*f)[i] = f_data + i * n;
  }
}

void deallocate_grids(float ***uo, float ***u, float ***f) {
  // de-allocate the data
  free((*u)[0]);
  free((*uo)[0]);
  free((*f)[0]);

  // de-allocate the rows pointers
  free(*u);
  free(*uo);
  free(*f);

  *u = NULL;
  *uo = NULL;
  *f = NULL;
}

void swap_grids(float ***uo, float ***u) {
  float **tmp = *uo;
  *uo = *u;
  *u = tmp;
}

void initialize_grids(int m, int n, int offset_m, int offset_n, float **uo, float **u, float **f, float h) {
  int prank, psize;
  MPI_Comm_rank(MPI_COMM_WORLD, &prank);
  MPI_Comm_size(MPI_COMM_WORLD, &psize);

  int i_start = (prank == 0 ? 0 : 1);
  int i_end   = (prank == psize - 1 ? m : m - 1);
  
  for (int i = i_start; i < i_end; i++) {
    for (int j = 0; j < n; j++) {
      u[i][j] = 0;
      uo[i][j] = 0;
      f[i][j] = -2. * 100. * M_PI * M_PI * sin(10. * M_PI * (i + offset_m - i_start) * h) *
          sin(10. * M_PI * (j + offset_n) * h);
    }
  }
}

static inline float compute_row(int i, int n, float **uo, float **u, float **f,
                                float h) {
  int j;
  float l2 = 0.;
  for (j = 1; j < n - 1; j++) {
    // computation of the new step
    u[i][j] = 0.25 * (uo[i - 1][j] + uo[i + 1][j] + uo[i][j - 1] +
                      uo[i][j + 1] - f[i][j] * h * h);

    // L2 norm
    l2 += (uo[i][j] - u[i][j]) * (uo[i][j] - u[i][j]);
  }
  return l2;
}

float compute_step(int m, int n, float **uo, float **u, float **f, float h) {
  float l2 = 0.;
    
  int prank, psize;
  MPI_Comm_rank(MPI_COMM_WORLD, &prank);
  MPI_Comm_size(MPI_COMM_WORLD, &psize);

  int north_prank = (prank == 0 ? MPI_PROC_NULL : prank - 1);
  int south_prank = (prank == (psize - 1) ? MPI_PROC_NULL : prank + 1);

  
  // Taking care of communications going up (so receiving from bottom)
  MPI_Sendrecv(uo[1], n, MPI_FLOAT, north_prank, 0,
               uo[m - 1], n, MPI_FLOAT, south_prank, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Taking care of communications going down (so receiving from top)
  MPI_Sendrecv(uo[m - 2], n, MPI_FLOAT, south_prank, 0,
               uo[0], n, MPI_FLOAT, north_prank, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  
  for (int i = 1; i < m - 1; i++) {
    l2 += compute_row(i, n, uo, u, f, h);
  }

  MPI_Allreduce(MPI_IN_PLACE, &l2, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  
  return l2;
}

float simulate(int m, int n, float **uo, float **u, float **f, float h, float epsilon, int * k) {
  float l2 = 0.;
  *k = 0;
  do {
    l2 = compute_step(m, n, uo, u, f, h);

    // copy new grid in old grid
    swap_grids(&uo, &u);

    ++(*k);
  } while (l2 > epsilon);

  // write the new old grid to disk
  write_to_file(m, n, uo, (*k), -1., 1.);

  return l2;
}
