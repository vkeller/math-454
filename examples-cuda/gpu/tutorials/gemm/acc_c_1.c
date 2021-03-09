#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#define M 1024 
#define N 2048 
#define P 4096 

/* muliply a M lines, P columns matrix by a P lines, N columns matrix 
   to produce a M lines, N columns matrix */ 
void 
MatrixMultiplication0(restrict float a[M][N], restrict float b[M][P], restrict float c[P][N]) 
{ 
    int i, j, k ; 

    for (i=0; i<M; i++) { 
        for (j=0; j<N; j++) { 
            for (k=0; k<P; k++) 
                a[i][j] += b[i][k]*c[k][j] ; 
        } 
    } 
} 

/* muliply a M lines, P columns matrix by a P lines, N columns matrix 
   to produce a M lines, N columns matrix */ 
void 
MatrixMultiplication1(restrict float a[M][N], restrict float b[M][P], restrict float c[P][N]) 
{ 
    int i, j, k ; 

#pragma acc region for parallel, vector(8) 
    for (i=0; i<M; i++) { 
#pragma acc for parallel, vector(8) 
        for (j=0; j<N; j++) { 
#pragma acc for seq 
            for (k=0; k<P; k++) 
                a[i][j] += b[i][k]*c[k][j] ; 
        } 
    } 
} 

void 
MatrixMultiplication2(restrict float a[M][N], restrict float b[M][P], restrict float c[P][N]) 
{ 
    int i, j, k ; 

#pragma acc region for parallel, vector(8) 
    for (i=0; i<M; i++){ 
#pragma acc for parallel, vector(8) 
        for (j=0; j<N; j++) { 
            float sum = 0.0 ; 
#pragma acc for seq 
            for (k=0; k<P; k++) 
                sum += b[i][k]*c[k][j] ; 
            a[i][j] = sum ; 
        } 
    } 
} 

void 
MatrixMultiplication3(float * restrict a, float * restrict b, float * restrict c, int m, int n, int p) 
{ 
    int i, j, k ; 

#pragma acc data region copyout(a[0:(m*n)-1]), copyin(b[0:(m*p)-1],c[0:(p*n)-1]) 
{ 
#pragma acc region for parallel, vector(8) 
    for (i=0; i<m; i++){ 
#pragma acc for parallel, vector (8) 
        for (j=0; j<n; j++) { 
#pragma acc for seq 
            for (k=0; k<p; k++) 
                a[i*n+j] += b[i*p+k]*c[k*n+j] ; 
   } 
    } 
} 
} 

void 
MatrixMultiplication4(float * restrict a,float * restrict b, float * restrict c, int m, int n, int p) 
{ 
    int i, j, k ; 

#pragma acc data region copyout(a[0:(m*n)-1]), copyin(b[0:(m*p)-1],c[0:(p*n)-1]) 
{ 
#pragma acc region for parallel, vector(8) 
    for (i=0; i<m; i++){ 
#pragma acc for parallel, vector (8) 
        for (j=0; j<n; j++) 
        { 
            float sum = 0.0 ; 
#pragma acc for seq 
            for (k=0; k<p; k++) 
                sum += b[i*p+k]*c[k*n+j] ; 
            a[i*n+j] = sum ; 
        } 
    } 
} 
} 

void main() 
{ 
    float *a, *b, *c; 
    int i; 
    float *a_1;
    struct timeval tv1, tv2;

    a = (float *) malloc(M*N*4); 
    a_1 = (float *) malloc(M*N*4); 
    b = (float *) malloc(M*P*4); 
    c = (float *) malloc(P*N*4); 

    for (i = 0; i <  M*N; i++) { 
      a[i] = (float) (i%3) -2.; 
      a_1[i] = a[i]; 
    } 
    for (i = 0; i <  M*P; i++) { 
      b[i] = (float) (i%5) -4.; 
    } 
    for (i = 0; i <  P*N; i++) { 
      c[i] = (float) (i%7)-5.; 
    } 

    printf("\tTesting simple C implementation of matrix multiplication.\n");
    gettimeofday(&tv1, NULL);
    MatrixMultiplication0((float **)a,(float **)b,(float **)c); 
    gettimeofday(&tv2, NULL);
    printf("\t\tdone... %f\n", a[10]);
    printf("\t\tExecution time (in millisec): %.2f\n",
	 (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
	 (double)(tv2.tv_sec -tv1.tv_sec )*1000);


    printf("\tTesting acc #1 of matrix multiplication.\n");
    gettimeofday(&tv1, NULL);
    MatrixMultiplication1((float **)a,(float **)b,(float **)c); 
    gettimeofday(&tv2, NULL);
    printf("\t\tdone... %f\n", a[10]);
    printf("\t\tExecution time (in millisec): %.2f\n",
	 (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
	 (double)(tv2.tv_sec -tv1.tv_sec )*1000);



      //MatrixMultiplication2((float **)a,(float **)b,(float **)c); 
      //MatrixMultiplication3((float *)a,(float *)b,(float *)c); 
      //MatrixMultiplication4((float *)a,(float *)b,(float *)c); 

 
    free(a); 
    free(b); 
    free(c); 
} 
