// --------------------------------------------------------------
// File: ex10.c
// Description: Example with the sections directive
// Software stack: Exercises suite for V. Keller lecture
// Version: 1.0
// License: BSD
// Author: Vincent Keller (Vincent.Keller@epfl.ch), CADMOS, 2012
//
// COMPILATION AND LINKING
//
// WITHOUT MKL :  gcc -O3 -ftree-vectorize dgemm.c -lgsl -lgslcblas -lm -o dgemm
// WITH MKL : export MKLROOT=/software/intel/mkl; icc -DMKLINUSE -DMKL_ILP64 -openmp -I${MKLROOT}/include -O3 -xHost dgemm.c  -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -o dgemm
// --------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#if defined(MKLINUSE)
#include "mkl.h"
#else
#include <gsl/gsl_blas.h>
#endif


double verification(double ** array, int N){
        double ret;
        int i,j;
        for (i=0;i<N;i++){
                for (j=0;j<N;j++){
                        ret+=array[i][j];
                }
        }
        return ret;
}

double second()
{
        struct timeval tp;
        int i;
        i = gettimeofday(&tp,NULL);
        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}


int main( )
{
        int N, myrank,i,j,k,chunk, mysize;
        double ** A;
        double ** B;
        double ** C;
        double t1, t2;

        double *A_dgemm, *B_dgemm, *C_dgemm;
        double alpha, beta;

        mysize = 1;
        N=2000;

#if defined(MKLINUSE)
        A_dgemm = (double *)mkl_malloc( N*N*sizeof( double ), 64 );
        B_dgemm = (double *)mkl_malloc( N*N*sizeof( double ), 64 );
        C_dgemm = (double *)mkl_malloc( N*N*sizeof( double ), 64 );
#else
        A_dgemm = (double *)malloc( N*N*sizeof( double ));
        B_dgemm = (double *)malloc( N*N*sizeof( double ));
        C_dgemm = (double *)malloc( N*N*sizeof( double ));
#endif

        for (i = 0; i < (N*N); i++) {
                A_dgemm[i] = 1.0;
        }

        for (i = 0; i < (N*N); i++) {
                B_dgemm[i] = 2.0;
        }

        for (i = 0; i < (N*N); i++) {
                C_dgemm[i] = 0.0;
        }

        alpha= 1.0;
        beta = 0.0;

        A = (double**)malloc(N*sizeof(double*));
        B = (double**)malloc(N*sizeof(double*));
        C = (double**)malloc(N*sizeof(double*));

        for (i=0;i<N;i++){
                A[i]=(double*)malloc(N*sizeof(double));
                B[i]=(double*)malloc(N*sizeof(double));
                C[i]=(double*)malloc(N*sizeof(double));
        }

        for (i=0;i<N;i++){
                for (j=0;j<N;j++){
                        A[i][j] = 1.0;
                        B[i][j] = 2.0;
                        C[i][j] = 0.0;
                }
        }

// ============================================================================================================
        t1 = second();
        for (i=0;i<N;i++){
                for (j=0;j<N;j++){
                        for (k=0;k<N;k++){
                                C[i][j]=C[i][j] + A[i][k]*B[k][j];
                        }
                }
        }
        t2 = second();

        printf("[IJK]   Compute time [s]   : %6.3f \n", (t2-t1));
        printf("[IJK]   Performance  [GF/s]: %6.3f \n", (((double)N*(double)N*(double)N*2.0)/((t2-t1)*1.e9)));
        printf("[IJK]   Verification       : %6.3f \n\n", verification(C,N));


// ============================================================================================================
        for (i=0;i<N;i++){
                for (j=0;j<N;j++){
                        A[i][j] = 1.0;
                        B[i][j] = 2.0;
                        C[i][j] = 0.0;
                }
        }
        t1 = second();
        for (i=0;i<N;i++){
                for (k=0;k<N;k++){
                        for (j=0;j<N;j++){
                                C[i][j]=C[i][j] + A[i][k]*B[k][j];
                        }
                }
        }
        t2 = second();

        printf("[IKJ]   Compute time [s]   : %6.3f \n", (t2-t1));
        printf("[IKJ]   Performance  [GF/s]: %6.3f \n", (((double)N*(double)N*(double)N*2.0)/((t2-t1)*1.e9)));
        printf("[IKJ]   Verification       : %6.3f \n\n", verification(C,N));


// ============================================================================================================
        for (i=0;i<N;i++){
                for (j=0;j<N;j++){
                        A[i][j] = 1.0;
                        B[i][j] = 2.0;
                        C[i][j] = 0.0;
                }
        }
        t1 = second();
        for (j=0;j<N;j++){
                for (i=0;i<N;i++){
                        for (k=0;k<N;k++){
                                C[i][j]=C[i][j] + A[i][k]*B[k][j];
                        }
                }
        }
        t2 = second();

        printf("[JIK]   Compute time [s]   : %6.3f \n", (t2-t1));
        printf("[JIK]   Performance  [GF/s]: %6.3f \n", (((double)N*(double)N*(double)N*2.0)/((t2-t1)*1.e9)));
        printf("[JIK]   Verification       : %6.3f \n\n", verification(C,N));


// ============================================================================================================
        for (i=0;i<N;i++){
                for (j=0;j<N;j++){
                        A[i][j] = 1.0;
                        B[i][j] = 2.0;
                        C[i][j] = 0.0;
                }
        }
        t1 = second();
        for (j=0;j<N;j++){
                for (k=0;k<N;k++){
                        for (i=0;i<N;i++){
                                C[i][j]=C[i][j] + A[i][k]*B[k][j];
                        }
                }
        }
        t2 = second();

        printf("[JKI]   Compute time [s]   : %6.3f \n", (t2-t1));
        printf("[JKI]   Performance  [GF/s]: %6.3f \n", (((double)N*(double)N*(double)N*2.0)/((t2-t1)*1.e9)));
        printf("[JKI]   Verification       : %6.3f \n\n", verification(C,N));


// ============================================================================================================
        for (i=0;i<N;i++){
                for (j=0;j<N;j++){
                        A[i][j] = 1.0;
                        B[i][j] = 2.0;
                        C[i][j] = 0.0;
                }
        }
        t1 = second();
        for (k=0;k<N;k++){
                for (i=0;i<N;i++){
                        for (j=0;j<N;j++){
                                C[i][j]=C[i][j] + A[i][k]*B[k][j];
                        }
                }
        }
        t2 = second();

        printf("[KIJ]   Compute time [s]   : %6.3f \n", (t2-t1));
        printf("[KIJ]   Performance  [GF/s]: %6.3f \n", (((double)N*(double)N*(double)N*2.0)/((t2-t1)*1.e9)));
        printf("[KIJ]   Verification       : %6.3f \n\n", verification(C,N));


// ============================================================================================================
        for (i=0;i<N;i++){
                for (j=0;j<N;j++){
                        A[i][j] = 1.0;
                        B[i][j] = 2.0;
                        C[i][j] = 0.0;
                }
        }
        t1 = second();
        for (k=0;k<N;k++){
                for (j=0;j<N;j++){
                        for (i=0;i<N;i++){
                                C[i][j]=C[i][j] + A[i][k]*B[k][j];
                        }
                }
        }
        t2 = second();

        printf("[KJI]   Compute time [s]   : %6.3f \n", (t2-t1));
        printf("[KJI]   Performance  [GF/s]: %6.3f \n", (((double)N*(double)N*(double)N*2.0)/((t2-t1)*1.e9)));
        printf("[KJI]   Verification       : %6.3f \n\n", verification(C,N));



// ============================================================================================================
        t1 = second();
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N,N,N, alpha, A_dgemm, N, B_dgemm, N, beta, C_dgemm, N);
        t2 = second();

        printf("[DGEMM OPT]   Compute time [s]   : %6.3f \n", (t2-t1));
        printf("[DGEMM OPT]   Performance  [GF/s]: %6.3f \n", (((double)N*(double)N*(double)N*2.0)/((t2-t1)*1.e9)));
        printf("[DGEMM OPT]   Verification       : %6.3f \n\n", verification(C,N));

        for (i=0;i<N;i++){
                free(A[i]);
                free(B[i]);
                free(C[i]);
        }

        free(A);
        free(B);
        free(C);

        return 0;
}
