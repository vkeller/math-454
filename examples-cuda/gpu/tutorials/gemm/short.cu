/* GEMM is a General Matrix Multiply - a subroutine in the Basic Linear Algebra Subprograms library*/

/* Includes, system */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

/* Includes, cuda */
//#include <cuda_runtime.h>
//#include <cublas_v2.h>
//#include <helper_cuda.h>

#define BLOCK_SIZE 16

/* ======================================================= */
/* CUDA implementation of dGEMM without using shared memory
/* ======================================================= */
__global__ void cuda_dgemm(int n, 
			   double alpha, 
			   const double *A, 
			   const double *B,
			   double beta, 
			   double *C) {

  int row = blockDim.y * blockIdx.y + threadIdx.y;
  int col = blockDim.x * blockIdx.x + threadIdx.x;
  
  printf("row = %d col = %d  n= %d\n", row, col, n);
  if (row > n || col > n) return;
  
  double prod = 0;
  for (int k = 0; k < n; ++k){
    prod += A[row * n + k] * B[k * n + col];
  }
  printf("alpha = %f, prod = %f\n", alpha, prod);
  C[row*n + col] = alpha * prod + beta * C[row*n+col]; 
}

/* ==== */
/* Main */
/* ==== */
int main(int argc, char **argv)
{
  //cublasStatus_t status;
  double *h_A, *h_B, *h_C, *h_C_blas, *h_C_simple, *h_C_0;
  double *d_A = 0; 
  double *d_B = 0;
  double *d_C = 0;
  double alpha = 1.0f;
  double beta = 0.0f;
  int n2, N;
  int i;
  double error_norm1, error_norm2;
  double ref_norm;
  double diff1, diff2;
  //cublasHandle_t handle;
  struct timeval tv1, tv2;


  /* get the size of the matrix from the command line */
  if (argc <2 ) N= 275;
  else N = atoi(argv[1]);
  n2 = N * N;

  /* Allocate host memory for the matrices */
  h_A = (double *)malloc(n2 * sizeof(double) );
  h_B = (double *)malloc(n2 * sizeof(double) );
  h_C = (double *)malloc(n2 * sizeof(double) );
  h_C_blas = (double *)malloc(n2 * sizeof(double) );
  h_C_simple = (double *)malloc(n2 * sizeof(double) );
  h_C_0 = (double *)malloc(n2 * sizeof(double) );

  /* Fill the matrices with test data */
  for (i = 0; i < n2; i++){
    h_A[i] = i; //rand() / (double)RAND_MAX;
    h_B[i] = i; //rand() / (double)RAND_MAX;
    h_C[i] = i; //rand() / (double)RAND_MAX;
    h_C_blas[i] = h_C[i];
    h_C_simple[i] = h_C[i];
    h_C_0[i] = h_C[i];
    printf("%f %f \n",h_A[i], h_B[i]);
  }

  /* ============ CUDA implementation without shared memory =============== */
  printf("\tTesting CUDA dgemm function without using Shared memory.\n");
  gettimeofday(&tv1, NULL);

  /* Allocate device memory for the matrices */
  cudaMalloc((void **)&d_A, n2 * sizeof(d_A[0]));
  cudaMalloc((void **)&d_B, n2 * sizeof(d_B[0]));
  cudaMalloc((void **)&d_C, n2 * sizeof(d_C[0]));

  /* copy A and B matrices to gpu */
  cudaMemcpy(d_A, h_A,n2*sizeof(d_A[0]), cudaMemcpyHostToDevice);
  cudaMemcpy(d_B, h_B,n2*sizeof(d_B[0]), cudaMemcpyHostToDevice);
  cudaMemcpy(d_C, h_C_0,n2*sizeof(d_C[0]), cudaMemcpyHostToDevice);

  /* Kernel */
  dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
  dim3 dimGrid(N/dimBlock.x+1, N/dimBlock.y+1);
  printf(" beta=%f\n",beta);
  cuda_dgemm<<<dimGrid, dimBlock>>>(N, alpha, d_A, d_B, beta, d_C);
  /* wait until all threads finish their job */
  cudaDeviceSynchronize();


  /* Read the result back */
  cudaMemcpy(h_C, d_C,n2*sizeof(d_C[0]), cudaMemcpyDeviceToHost);


  gettimeofday(&tv2, NULL);
  printf("\t\tdone...\n");
  printf("\t\tExecution time (in millisec): %.2f\n",
	 (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
	 (double)(tv2.tv_sec -tv1.tv_sec )*1000);

  printf("\n\tChecking results.\n");
  /* Check result against reference */
  error_norm1 = 0;
  ref_norm = 0;
  for (i = 0; i < n2; ++i){
    if (i<10)printf("%f\n",h_C[i]);

  }
  /* free cuda memory */
  cudaFree(d_A);
  cudaFree(d_B);
  cudaFree(d_C);

  /* Memory clean up */
  free(h_A);
  free(h_B);
  free(h_C);
  free(h_C_simple);
  free(h_C_blas);

  return(0);
}
