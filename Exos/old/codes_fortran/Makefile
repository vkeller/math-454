F90		= ifort
CC		= icc
MPIF90		= mpiifort
FFLAGS		= -O3 -xHost -c -diag-enable=vec -g
FFLAGS_OMP 	= $(FFLAGS) -openmp -openmp-report=2
FFLAGS_OFFLOAD 	= $(FFLAGS) -offload -offload-attribute-target=mic OFFLOAD_REPORT=2
FFLAGS_MIC 	= -mmic -O3
LDFLAGS		= write_to_file_binary1D.o poisson_module.o -g -lm
LDFLAGS_OMP	= $(LDFLAGS) -openmp 

EXEC_HOST_1D	= poisson_host_f_1D
EXEC_HOST_2D	= poisson_host_f_2D
EXEC_OMP	= poisson_omp_f
EXEC_MPI_2D	= poisson_mpi_f_2D
EXEC_MPI_1D	= poisson_mpi_f_1D
EXEC_MPI_P	= poisson_mpi_persistent_f
EXEC_MPI_P_IO	= poisson_mpi_persistent_mpiio_f
EXEC_MPIIO	= poisson_mpiio_f
EXEC_OFFLOAD	= poisson_offload_f
EXEC_MIC	= poisson_mic_f
EXEC_CUDA	= poisson_cuda_f




#all: poisson_host poisson_omp poisson_mpi poisson_mpiio poisson_offload poisson_mic
all: poisson_host_2D poisson_host_1D poisson_mpi_1D poisson_mpi_2D


poisson_host_1D: 
	$(CC) -c -g write_to_file_binary1D.c
	$(F90) $(FFLAGS) poisson_module.f90
	$(F90) $(FFLAGS) poisson1D.f90
	$(F90) $(LDFLAGS) poisson1D.o -o $(EXEC_HOST_1D)

poisson_host_2D: 
	$(CC) -c -g write_to_file_binary1D.c
	$(F90) $(FFLAGS) poisson_module.f90
	$(F90) $(FFLAGS) poisson2D.f90
	$(F90) $(LDFLAGS) poisson2D.o -o $(EXEC_HOST_2D)

poisson_mpi_1D: 
	$(CC) -c -g write_to_file_binary1D.c
	$(F90) $(FFLAGS) poisson_module.f90
	$(MPIF90) $(FFLAGS) poisson_mpi1D.f90
	$(MPIF90) $(LDFLAGS) poisson_mpi1D.o -o $(EXEC_MPI_1D)


poisson_mpi_2D: 
	$(CC) -c -g write_to_file_binary1D.c
	$(F90) $(FFLAGS) poisson_module.f90
	$(MPIF90) $(FFLAGS) poisson_mpi2D.f90
	$(MPIF90) $(LDFLAGS) poisson_mpi2D.o -o $(EXEC_MPI_2D)


poisson_omp: 
	$(F90) $(FFLAGS_OMP) poisson_omp.f90
	$(F90) $(LDFLAGS_OMP) poisson_omp.o -o $(EXEC_OMP)

poisson_mpi_persistent: 
	$(MPIF90) $(FFLAGS) poisson_mpi_persistent.f90
	$(MPIF90) $(LDFLAGS) poisson_mpi_persistent.o -o $(EXEC_MPI_P)

poisson_mpi_persistent_mpiio: 
	$(MPIF90) $(FFLAGS) poisson_mpi_persistent_mpiio.f90
	$(MPIF90) $(LDFLAGS) poisson_mpi_persistent_mpiio.o -o $(EXEC_MPI_P_IO)

poisson_mpiio: 
	$(MPIF90) $(FFLAGS) poisson_mpiio.f90
	$(MPIF90) $(LDFLAGS) poisson_mpiio.o -o $(EXEC_MPIIO)

poisson_offload: 
	$(F90) $(FFLAGS_OFFLOAD) poisson_offload.f90
	$(F90) $(LDFLAGS) poisson_offload.o -o $(EXEC_OFFLOAD)

poisson_mic: 
	$(F90) $(FFLAGS_MIC) poisson_mic.f90
	$(F90) $(LDFLAGS) poisson_mic.o -o $(EXEC_MIC)

poisson_cuda: 
	$(F90) $(FFLAGS_MIC) poisson_mic.f90
	$(F90) $(LDFLAGS) poisson_mic.o -o $(EXEC_MIC)

clean:
	rm -f *.bmp *.o *.mod $(EXEC_HOST_1D) $(EXEC_HOST_2D) $(EXEC_MPI_1D) $(EXEC_MPI_2D) $(EXEC_MPIIO) $(EXEC_OFFLOAD) $(EXEC_MIC) $(EXEC_OMP) $(EXEC_MPI_P) $(EXEC_MPI_P_IO)  
