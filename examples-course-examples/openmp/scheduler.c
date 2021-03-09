#include "timing.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

double verification(double ** matrix, int N){
  double ret = 0.;
  int i,j;

  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      ret+=matrix[i][j];
    }
  }
  
  return ret;
}

int main() {
  int N, i, j, k;
  int mysize, chunk;
  double ** A, ** B, ** C;
  double * A_data, * B_data, * C_data;
  double t1, t2;

  N=2000;
  
  A_data = (double *)malloc(N*N * sizeof(double));
  B_data = (double *)malloc(N*N * sizeof(double));
  C_data = (double *)malloc(N*N * sizeof(double));

  A = (double **)malloc(N * sizeof(double *));
  B = (double **)malloc(N * sizeof(double *));
  C = (double **)malloc(N * sizeof(double *));

  if(A_data == NULL) {
    printf("Cannot allocate A\n");
    return 1;
  }

  if(B_data == NULL) {
    printf("Cannot allocate B\n");
    return 1;
  }

  if(C_data == NULL) {
    printf("Cannot allocate C\n");
    return 1;
  }

  for (i = 0; i < (N*N); i++) {
    A_data[i] = 1.0;
    B_data[i] = 2.0;
    C_data[i] = 0.0;
  }
  
  for (i = 0; i < N; i++){
    A[i] = A_data + i * N;
    B[i] = B_data + i * N;
    C[i] = C_data + i * N;
  }

  t1 = second();
#pragma omp parallel shared(A, B, C) private(i, j, k)
  {
    mysize = omp_get_num_threads();
    chunk = (N/mysize);

#pragma omp for schedule(static, chunk)
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
	for (k = 0; k < N; k++) {
	  C[i][j] = C[i][j] + A[i][k] * B[k][j];
	}
      }
    }
  }
  
  t2 = second();
  
  printf("[IJK]   Nb cores [-]       : %i \n", (mysize));
  printf("[IJK]   Compute time [s]   : %6.3f \n", (t2-t1));
  printf("[IJK]   Performance  [GF/s]: %6.3f \n", ((pow((double)N, 3.)*2.0)/((t2-t1)*1.e9)));
  printf("[IJK]   Verification       : %6.3f \n\n", verification(C,N));
  
  free(A);
  free(B);
  free(C);

  free(A_data);
  free(B_data);
  free(C_data);

  return 0;
}
