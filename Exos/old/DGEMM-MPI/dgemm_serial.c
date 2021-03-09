#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#define n 1024

#define ID(i, j) ((i) * n + (j))

int main(int argc, char *argv[]) {
	int     i, j, k;
	float **a;
	float **b;
	float **c;

	double start, end;
	double flops,gflops,elapsedtime;
	struct timeval tv1,tv2;

  // Allocation of arrays of pointers to the number of rows per proc
	a  = (float**) malloc(n * sizeof(float*));
	b = (float**) malloc(n * sizeof(float*));
	c  = (float**) malloc(n * sizeof(float*));

  // allocation of the rows to the size of full lines
	for(i=0; i < n; i++) {
		a [i] = (float*) malloc(n*sizeof(float));
		b [i] = (float*) malloc(n*sizeof(float));
		c [i] = (float*) malloc(n*sizeof(float));
	}

	for (i=0;i<n;i++){
		for (j=0;j<n;j++){
			a[i][j] = 1.; b[i][j] = 1.; c[i][j] = 1.;
		}
	}


	gettimeofday(&tv1, (struct timezone*)0);
	for (i=0;i<n;i++){
		for (k=0;k<n;k++){
			for (j=0;j<n;j++){
				c[i][j]=c[i][j]+a[i][k]*b[k][j];
			}
		}
	}
	gettimeofday(&tv2, (struct timezone*)0);
	elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
	flops=(double)n*(double)n*(double)n*2.0;
	gflops=flops/(elapsedtime*pow((double)10,9));
	printf("Gflops = %E \t elaspedtime = %E\t size = %d\n",gflops,elapsedtime,n);

  // deallocate the pointers
	for(i=0; i < n; i++) {
		free(a[i]);free(b[i]);free(c[i]);
	}  
	free(a);free(b);free(c);

  return 0;
}

