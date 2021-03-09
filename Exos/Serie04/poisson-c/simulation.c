/* -------------------------------------------------------------------------- */
#include "io_binary.h"
/* -------------------------------------------------------------------------- */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/* -------------------------------------------------------------------------- */

void allocate_grids(int n, float ***u, float ***uo, float ***f) {
  int i;

  // Allocation of the data storage
  float *u_data = (float *)malloc(n * n * sizeof(float));
  float *uo_data = (float *)malloc(n * n * sizeof(float));
  float *f_data = (float *)malloc(n * n * sizeof(float));

  // Allocation of arrays of pointers for rows
  *u = (float **)malloc(n * sizeof(float *));
  *uo = (float **)malloc(n * sizeof(float *));
  *f = (float **)malloc(n * sizeof(float *));

  // set the row pointers in the memory
  for (i = 0; i < n; i++) {
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

void initialize_grids(int n, float **uo, float **u, float **f, float h) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      u[i][j] = 0;
      uo[i][j] = 0;
      f[i][j] = -2. * 100. * M_PI * M_PI * sin(10. * M_PI * i * h) *
                sin(10. * M_PI * j * h);
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

float compute_step(int n, float **uo, float **u, float **f, float h) {
  float l2 = 0.;
  int i;

  for (i = 1; i < n - 1; i++) {
    l2 += compute_row(i, n, uo, u, f, h);
  }

  return l2;
}

float simulate(int n, float **uo, float **u, float **f, float h, float epsilon, int * k) {
  float l2 = 0.;
  *k = 0;
  do {
    l2 = compute_step(n, uo, u, f, h);

    // copy new grid in old grid
    swap_grids(&uo, &u);

    // write the new old grid to disk
    write_to_file(n, uo, (*k), -1., 1.);

    ++(*k);
  } while (l2 > epsilon);
  
  return l2;
}
