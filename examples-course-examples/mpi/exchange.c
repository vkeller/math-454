#include <mpi.h>
#include <string.h>
#include <stdio.h>

int main(int argc, char * argv[]) {
  int myrank;
  int buf1[100];
  int buf2[100];

  MPI_Status status;
  MPI_Request request;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if (myrank == 0) {
    MPI_Isend(buf1, 100, MPI_INT, 1, 0, MPI_COMM_WORLD, &request);
    MPI_Recv(buf2, 100, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else if (myrank == 1) {
    MPI_Isend(buf1, 100, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
    MPI_Recv(buf2, 100, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  MPI_Wait(&request, &status);
  memcpy(buf1, buf2, 100 * sizeof(int));

  MPI_Finalize();

  return 0;
}
