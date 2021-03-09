/* GEMM is a General Matrix Multiply - a subroutine in the Basic Linear Algebra Subprograms library*/

/* Includes, system */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>


/* ======================================================= */
/* Simple implementation of dgemm */
/* ======================================================= */
static void simple_dgemm(int n, double alpha, const double *A, const double *B,
                         double beta, double *C) {
  int i, j, k;
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j){
      double prod = 0;
      for (k = 0; k < n; ++k){
	prod += A[k * n + i] * B[j * n + k];
      }
      C[j * n + i] = alpha * prod + beta * C[j * n + i];
    }
  }
}

/* ======================================================= */
/* dgemm implementation with openACC acceleration*/
/* ======================================================= */
static void acc_dgemm(int n, double alpha, const double *A, const double *B,
                         double beta, double *C) {
  int i, j, k;

  #pragma acc parallel loop copyin(A[0:(n*n)], B[0:(n*n)]) copy(C[0:(n*n)])
  for (i = 0; i < n; ++i) {

    #pragma acc loop
    for (j = 0; j < n; ++j){

      double prod = 0;
      for (k = 0; k < n; ++k) prod += A[k * n + i] * B[j * n + k];
      
      C[j * n + i] = alpha * prod + beta * C[j * n + i];
    }
  }
}



/* ==== */
/* Main */
/* ==== */
int main(int argc, char **argv)
{
  double *h_A, *h_B, *h_C, *h_C_simple;
  double alpha = 1.0f;
  double beta = 0.0f;
  int n2, N;
  int i;
  double error_norm;
  double ref_norm;
  double diff;

  struct timeval tv1, tv2;


  /* get the size of the matrix from the command line */
  if (argc <2 ) N= 275;
  else N = atoi(argv[1]);
  n2 = N * N;

  printf("\nRunning dgemm test for %d by %d matricies.\n", N, N);
  
  /* Allocate host memory for the matrices */
  h_A = (double *)malloc(n2 * sizeof(double) );
  h_B = (double *)malloc(n2 * sizeof(double) );
  h_C = (double *)malloc(n2 * sizeof(double) );
  h_C_simple = (double *)malloc(n2 * sizeof(double) );

  /* Fill the matrices with test data */
  for (i = 0; i < n2; i++){
    h_A[i] = rand() / (double)RAND_MAX;
    h_B[i] = rand() / (double)RAND_MAX;
    h_C[i] = rand() / (double)RAND_MAX;
    h_C_simple[i] = h_C[i];
  }

  printf("\tTesting simple C implementation of dgemm function.\n");
  gettimeofday(&tv1, NULL);
  /* Performs operation using plain C code */
  simple_dgemm(N, alpha, h_A, h_B, beta, h_C_simple);
  gettimeofday(&tv2, NULL);
  printf("\t\tdone...\n");
  printf("\t\tExecution time (in millisec): %.2f\n",
	 (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
	 (double)(tv2.tv_sec -tv1.tv_sec )*1000);


  printf("\tTesting dgemm function using openACC directives.\n");
  gettimeofday(&tv1, NULL);
  acc_dgemm(N, alpha, h_A, h_B, beta, h_C);

  gettimeofday(&tv2, NULL);
  printf("\t\tdone...\n");
  printf("\t\tExecution time (in millisec): %.2f\n",
	 (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
	 (double)(tv2.tv_sec -tv1.tv_sec )*1000);

  printf("\n\tChecking results.\n");
  /* Check result against reference */
  error_norm = 0;
  ref_norm = 0;
  for (i = 0; i < n2; ++i){
    diff = h_C_simple[i] - h_C[i];
    error_norm += diff * diff;
  }

  error_norm = (double)sqrt((double)error_norm);

  if(error_norm  < 1e-6f)printf("\t\tPassed acc Dgemm vs. simple Dgemm comparison!\n");
  else printf("\t\tDid not pass acc Dgemm vs. simple Dgemm comparison!\n");
  
  /* Memory clean up */
  free(h_A);
  free(h_B);
  free(h_C);
  free(h_C_simple);

  return(0);
}
