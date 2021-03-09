#include <mpi.h>
#include <stdio.h>

int main(int argc, char * argv[]) {
  int rank;
  int buf[100];

  MPI_Request request;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
    MPI_Isend(buf, 100, MPI_INT, 1, 0, MPI_COMM_WORLD, &request);
  else if (rank == 1)
    MPI_Irecv(buf, 100, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);

  MPI_Wait(&request, &status);
  MPI_Finalize();
}
