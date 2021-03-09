
// cuda-kernel: add 2 numbers 
__global__ void addnums (double *pi, double c){

  *pi += c;
}

// cuda-kernel: add 2 vectors
__global__ void addvecs (double *v1, double *v2){

  int idx = threadIdx.x;
    v1[idx] += v2[idx];
}
