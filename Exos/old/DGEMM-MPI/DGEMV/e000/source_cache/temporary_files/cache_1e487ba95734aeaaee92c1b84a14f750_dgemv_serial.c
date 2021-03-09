#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#define n 8
#define ntime 2000

#define ID(i, j) ((i) * n + (j))

int main(int argc, char *argv[]) {
	int     i, j, k, ndiag, nn;
	float **A;
	float *Ua;
	float *U;

	double start, end;
	double flops,gflops,elapsedtime;
	struct timeval tv1,tv2;


	for (k=1;k<=ntime;k++){
		nn = k*n;
		ndiag = nn;
// Allocation of arrays of pointers to the number of rows per proc
		A  = (float**) malloc(nn * sizeof(float*));
		Ua = (float*) malloc(ndiag * sizeof(float*));
		U  = (float*) malloc(ndiag * sizeof(float*));

// allocation of the rows to the size of full lines
		for(i=0; i < nn; i++) {
			A [i] = (float*) malloc(nn*sizeof(float));
		}

		for (i=0;i<ndiag;i++){
			U[i] = 2.; Ua[i] = 0.;
			for (j=0;j<nn;j++){
				A[i][j] = 3.; 
			}
		}


		gettimeofday(&tv1, (struct timezone*)0);	
		for (i=0;i<nn;i++){
			for (j=0;j<ndiag;j++){
				Ua[i] = Ua[i] + A[i][j]*U[j];
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		elapsedtime=elapsedtime/(double)ntime;
		flops=(double)n*(double)n*(double)n*2.0;
		gflops = 2.*(double)n*(double)ndiag / elapsedtime / 1.0E9;
//	gflops=flops/(elapsedtime*pow((double)10,9));
//	printf("Gflops = %E \t elaspedtime = %E\t size = %d\n",gflops,elapsedtime,n);
//		printf("[k = %d] Gflops = %E \t elaspedtime = %E\t size = %d\n",k,gflops,elapsedtime,nn);
		printf("%d\t %E\n",nn,gflops);

// deallocate the pointers
		for(i=0; i < nn; i++) {
			free(A[i]);
		}  
		free(A);free(Ua);free(U);
	}

	return 0;
}

