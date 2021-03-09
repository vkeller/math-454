#include <mpi.h>
#include <stdio.h>

int main(int argc, char * argv[]) {
  int myrank, mysize;
  int buf[100];
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &mysize);

  if (mysize != 2)
    printf("Warning: this examples will most probably deadlock with a number "
           "of process different from 2.\n");

  if (myrank == 0) {
    MPI_Send(buf, 100, MPI_INT, 1, 0, MPI_COMM_WORLD);
  } else {
    MPI_Recv(buf, 100, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
  }

  MPI_Finalize();
}
