#include <stdio.h>
#include <stdlib.h>
#include <omp.h>


int main(int argc, char *argv[]) {
  int * vec;
  int global_sum, i;
  int size_vec = 10;
  
  vec = (int*) malloc (size_vec*sizeof(int));
  global_sum = 0;

  for (i = 0; i < size_vec; i++) {
    vec[i] = i;
  }

#pragma omp parallel for reduction(+:global_sum)
  for (i = 0; i < size_vec; i++) {
    global_sum += vec[i];
  }

  printf("sum = %i\n", global_sum);

  return 0;
}
