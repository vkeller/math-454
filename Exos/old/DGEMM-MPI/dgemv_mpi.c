#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <mpi.h>

#define n 64
#define ntime 300

#define ID(i, j) ((i) * n + (j))

int main(int argc, char *argv[]) {
	int     i, j, k;
	int	ndiag, nn, nn_loc;
	float *A, *A_loc;
	float *Ua;
	float *U;

	double start, end;
	double flops,gflops;
	struct timeval tv1,tv2;
	double t1,t2,t3,t4,elapsedtime_comp, elapsedtime_comm, elapsedtime, total_time;

	int myrank, mysize;


	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &mysize);

	if (myrank==0) printf("MXV \nNumber of Processes : %d \n\nSize\t GFlops  \t time_comp\t time_comm\t (t_comp/t_comm)\n=====================================================================\n",mysize);

	total_time = 0.;

	for (k=1;k<=ntime;k++){
		nn = k*n;
		ndiag = nn;
		nn_loc = (nn / mysize) + (myrank < nn % mysize ? 1 : 0);
// Allocation of arrays of pointers to the number of rows per proc
// Size of A is A[nn:nn_loc]
/*
		A_loc  = (float*) aligned_alloc(sizeof(float*), nn_loc*ndiag);
		Ua = (float*) aligned_alloc(sizeof(float*), ndiag);
		U  = (float*) aligned_alloc(sizeof(float*), ndiag);
*/
		A_loc  = (float*) malloc(sizeof(float*)*nn_loc*ndiag);
		Ua = (float*) malloc(sizeof(float)*nn_loc);
		U  = (float*) malloc(sizeof(float)*ndiag);

// allocation of the rows to the size of full lines

		if (myrank==0) {
//			A  = (float*) aligned_alloc(sizeof(float*), nn_loc*ndiag);
			A  = (float*) malloc(sizeof(float)*ndiag*ndiag);
			for (j=0;j<ndiag;j++){
				U[i] = 2.;
				for (i=0;i<ndiag;i++){
					A[i+j*nn] = 3.; 
				}
			}
		}

		for (i=0;i<nn_loc;i++){
			Ua[i] = 0.;
		}

		t1 = MPI_Wtime();
		MPI_Bcast(U, ndiag, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Scatter(A,(nn_loc*ndiag),MPI_FLOAT,A_loc,(nn_loc*ndiag),MPI_FLOAT,0,MPI_COMM_WORLD);
		t2 = MPI_Wtime();		
		for (j=0;j<ndiag;j++){
                        for (i=0;i<nn_loc;i++){
                                Ua[i] = Ua[i] + A_loc[i+j*nn_loc]*U[i];
                        }
                }

		t3 = MPI_Wtime();
		MPI_Gather(Ua,nn_loc,MPI_FLOAT,U,nn_loc,MPI_FLOAT,0,MPI_COMM_WORLD);
		t4 = MPI_Wtime();
		elapsedtime = (t4-t1);
		elapsedtime_comm = (t2-t1) + (t4-t3);
		elapsedtime_comp = (t3-t2);	
		total_time += elapsedtime_comp + elapsedtime_comm;
		flops=(double)n*(double)n*(double)n*2.0;
		gflops = 2.*(double)n*(double)ndiag / elapsedtime / 1.0E9;
//		if (myrank==0) printf("%d\t %E Computation time = %E \t Communication time = %E \t Ratio(t_comp/t_comm) = %E\n",nn,gflops,elapsedtime_comp,elapsedtime_comm,(elapsedtime_comp/elapsedtime_comm));
		if (myrank==0) printf("%d\t %E\t %E\t %E\t %E\n",nn,gflops,elapsedtime_comp,elapsedtime_comm,(elapsedtime_comp/elapsedtime_comm));

// deallocate the pointers
		free(A_loc);free(Ua);free(U);
		if (myrank==0) free(A);
	}
	if (myrank==0) printf("\nTotal time : %E seconds\n",total_time);

	MPI_Finalize();

	return 0;
}

