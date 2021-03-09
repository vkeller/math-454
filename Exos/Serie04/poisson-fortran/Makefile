F90		= gfortran
CC		= gcc
MPIF90		= mpif90
#FFLAGS		= -O3 -xHost -c -diag-enable=vec -g
FFLAGS		= -O3 -ftree-vectorize -g -c -ffree-line-length-none
LDFLAGS		= write_to_file_binary1D.o poisson_module.o -g -lm

EXEC_HOST_1D	= poisson_host_f_1D

all: poisson_host_1D


poisson_host_1D: 
	$(CC) -c -g write_to_file_binary1D.c
	$(F90) $(FFLAGS) poisson_module.f90
	$(F90) $(FFLAGS) poisson1D.f90
	$(F90) $(LDFLAGS) poisson1D.o -o $(EXEC_HOST_1D)

clean:
	rm -f *.bmp *.o *.mod $(EXEC_HOST_1D)
