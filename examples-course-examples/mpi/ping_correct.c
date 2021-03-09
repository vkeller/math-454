#include <mpi.h>
#include <stdio.h>

int main(int argc, char * argv[]) {
  int myrank, mysize;
  int buf[100];

  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &mysize);

  if (myrank == 0) {
    MPI_Send(buf, 100, MPI_INT, 1, 0, MPI_COMM_WORLD);
  } else if (myrank == 1) {
    MPI_Recv(buf, 100, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
  }

  MPI_Finalize();
}
