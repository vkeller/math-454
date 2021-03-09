#include <stdio.h>
#include <stdlib.h>

#define up 10

int main() {
  int i, n;
  n = 0;

  for (i = 0; i < up; i++){
    n = n + 1;
  }

  printf("n = %i\n", n);

  return 0;
}
