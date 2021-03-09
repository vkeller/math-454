CC		= icc
MPICC		= mpiicc
CFLAGS		= -O3 -Wall -xHost -c -diag-enable=vec -g
CFLAGS_OMP 	= $(CFLAGS) -openmp -openmp-report=2
CFLAGS_OFFLOAD 	= $(CFLAGS) -offload -offload-attribute-target=mic OFFLOAD_REPORT=2
CFLAGS_MIC 	= -mmic -O3
LDFLAGS		= -g
LDFLAGS_OMP	= -openmp

EXEC_HOST	= poisson_host
EXEC_OMP	= poisson_omp
EXEC_MPI	= poisson_mpi
EXEC_MPI_P	= poisson_mpi_persistent
EXEC_MPI_P_IO	= poisson_mpi_persistent_mpiio
EXEC_MPIIO	= poisson_mpiio
EXEC_OFFLOAD	= poisson_offload
EXEC_MIC	= poisson_mic



#all: poisson_host poisson_omp poisson_mpi poisson_mpiio poisson_offload poisson_mic
all: poisson_host poisson_omp poisson_mpi poisson_mpiio


poisson_host: 
	$(CC) $(CFLAGS) poisson1D.c
	$(CC) $(LDFLAGS) poisson1D.o -o $(EXEC_HOST)

poisson_omp: 
	$(CC) $(CFLAGS_OMP) poisson_omp.c
	$(CC) $(LDFLAGS_OMP) poisson_omp.o -o $(EXEC_OMP)

poisson_mpi: 
	$(MPICC) $(CFLAGS) poisson_mpi.c
	$(MPICC) $(LDFLAGS) poisson_mpi.o -o $(EXEC_MPI)

poisson_mpi_persistent: 
	$(MPICC) $(CFLAGS) poisson_mpi_persistent.c
	$(MPICC) $(LDFLAGS) poisson_mpi_persistent.o -o $(EXEC_MPI_P)

poisson_mpi_persistent_mpiio: 
	$(MPICC) $(CFLAGS) poisson_mpi_persistent_mpiio.c
	$(MPICC) $(LDFLAGS) poisson_mpi_persistent_mpiio.o -o $(EXEC_MPI_P_IO)

poisson_mpiio: 
	$(MPICC) $(CFLAGS) poisson_mpiio.c
	$(MPICC) $(LDFLAGS) poisson_mpiio.o -o $(EXEC_MPIIO)

poisson_offload: 
	$(CC) $(CFLAGS_OFFLOAD) poisson_offload.c
	$(CC) $(LDFLAGS) poisson_offload.o -o $(EXEC_OFFLOAD)

poisson_mic: 
	$(CC) $(CFLAGS_MIC) poisson_mic.c
	$(CC) $(LDFLAGS) poisson_mic.o -o $(EXEC_MIC)

clean:
	rm -f *.o $(EXEC_HOST) $(EXEC_MPI) $(EXEC_MPIIO) $(EXEC_OFFLOAD) $(EXEC_MIC) $(EXEC_OMP) $(EXEC_MPI_P) $(EXEC_MPI_P_IO)  
