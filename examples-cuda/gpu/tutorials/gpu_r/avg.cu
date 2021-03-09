#include <stdlib.h>
#include <stdio.h>


__global__ void avg_kernel(double *heights, 
			   int nx, int ny, 
			   double Radius,
			   double *output){

  int i,j, ind;
  int ix, iy;
  int ixmin, ixmax, iymin, iymax;
  double h;
  int N;
  double ave;


  i = blockIdx.x * blockDim.x + threadIdx.x;
  j = blockIdx.y * blockDim.y + threadIdx.y; 

  // check the array boundaries
  if (i >= nx || j >= ny) return;

  ind = j * nx + i;  // location in the array
  h = heights[ind];
  output [ind] = 0;


  ixmin = max( i - int(Radius) , 0); 
  ixmax = min( i + int(Radius), nx-1); 
  iymin = max( j - int(Radius) , 0); 
  iymax = min( j + int(Radius), ny-1);

  N=0;
  ave = 0; 
  for(ix = ixmin; ix <= ixmax; ix++){
    for(iy = iymin; iy <= iymax; iy++){
      if ((ix-i)*(ix-i) + (iy-j)*(iy-j) <= Radius*Radius ){
	N++;
	ave = ave + heights[iy*nx+ix];
      }
    }
  }
  if (N > 0) output [ind] = ave/N;
  return;

}


#define BLOCK_SIZE 16

extern "C"{
  void avg_cuda(double *heights,  
		int nx,
		int ny,
		double r,
		double *output){


    dim3 nThreads(BLOCK_SIZE,BLOCK_SIZE);
    dim3 nBlocks ( (nx-1)/BLOCK_SIZE + 1, (ny-1)/BLOCK_SIZE + 1);

    double *d_heights;
    double *d_output;


    // allocate memory on GPU
    cudaMalloc((void**) &d_heights, nx*ny * sizeof(double));
    cudaMalloc((void**) &d_output, nx*ny * sizeof(double));

    // copy input array:
    cudaMemcpy(d_heights, heights, nx*ny*sizeof(double),cudaMemcpyHostToDevice);

    // execute Kernel
    avg_kernel<<<nBlocks,nThreads>>>(d_heights, nx, ny, r, d_output);

    // copy output array back to the CPU
    cudaMemcpy(output, d_output, nx*ny*sizeof(double),cudaMemcpyDeviceToHost);

    // free the memory
    cudaFree(d_heights);
    cudaFree(d_output);

  }

}

