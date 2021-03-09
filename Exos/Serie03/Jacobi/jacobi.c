#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <sys/time.h>
#include <math.h>

#include "utils.h"

#define MIN(a,b)  (((a<b)?a:b))

#define DIM_N 1000
#define DIM_M 1000

#define NREP  1
//
void laplacian(double*, double*, int, int);

static double EPSI=1.0e-1;

int main(int argc, char** argv) 
{
	//
	int dim_n = DIM_N;
	int dim_m = DIM_M;
	printf("argc = %d\n", argc);
	//      
	if (argc == 2)
	{
		dim_n = dim_m = atoi(argv[1]);
	}
	//
	double* restrict storage1 = (double*)malloc(dim_n*dim_m*sizeof(double));
	double* restrict storage2 = (double*)malloc(dim_n*dim_m*sizeof(double));
	//
	printf("3x3 stencil...%dx%d\n\n", dim_m, dim_n);
	printf("array sizes = %f MB\n", (double) dim_n*dim_m*sizeof(double)/1024/1024);
	//
	double alloctime = -myseconds();
	init (storage1, dim_n, dim_m);
	init (storage2, dim_n, dim_m);
	alloctime += myseconds();
	printf("Allocation time = %f s.\n\n", alloctime);
	//
	int nrep    = NREP;
	int ii;
	//
	double * st_read  = storage2;
	double * st_write = storage1;
	//
	double norminf, norml2;
	double time = - myseconds();
	int count = 0;	
	//
	do 
	{
		count++;
		// swaping the solutions
		double *tmp = st_read;
		st_read     = st_write;
		st_write    = tmp;
		//
		double ftime = -myseconds();
		unsigned long long c0 = cycles();
		//
		for (ii = 0; ii < nrep; ++ii)
		{
			laplacian(st_read, st_write, dim_m, dim_n);
		}
		//
		unsigned long long c1 = cycles();
		unsigned long long cycles = (c1 - c0)/nrep;
		//
		ftime += myseconds();
		ftime /= nrep;
		//
		norminf  = maxNorm(st_write, st_read, dim_n*dim_m);
		norml2   = l2Norm (st_write, st_read, dim_n*dim_m);
		//
		double flops = (dim_m - 2)*(dim_n - 2)*4.;
		//
		printf("iter = %d, linf norm = %g l2 norm = %g, %d flops, %ld cycles, %f flops/cycle, %f s., %f Gflops/s\n", count, norminf, norml2, (int) flops, cycles, flops/cycles, ftime, (double) flops/ftime/1e9); 
	} while (count < 20);

	time += myseconds();

	printf("\n# iter= %d time= %g residual= %g\n", count, time, norml2); 

	free(storage1);
	free(storage2);

	return 0;
}
