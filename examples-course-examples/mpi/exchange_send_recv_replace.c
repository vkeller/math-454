#include <mpi.h>
#include <string.h>
#include <stdio.h>

int main(int argc, char * argv[]) {
  int myrank;
  int buf1[100];

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if (myrank == 0) {
    MPI_Sendrecv_replace(buf1, 100, MPI_INT, 1, 0, 1, 0, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
  } else if (myrank == 1) {
    MPI_Sendrecv_replace(buf1, 100, MPI_INT, 0, 0, 0, 0, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
  }

  MPI_Finalize();

  return 0;
}
