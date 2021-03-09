
// cuda-kernel: add 2 vectors
__global__ void addvecs (double *v1, double *v2){

  int idx = threadIdx.x;
    v1[idx] += v2[idx];
}
