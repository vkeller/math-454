GSL_LIB=-L${GSL_LIBRARY} -lgsl -lgslcblas
GSL_INC=-I${GSL_INCLUDE}
GCC_OPT=-O3 -ftree-vectorize

MKL_FLAGS=-DMKLINUSE -DMKL_ILP64 -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread
ICC_OPT=-O3 -xHost

dgemm_gcc: dgemm.c
	gcc $(GCC_OPT) dgemm.c  $(GSL_LIB) $(GSL_INC) -lm -o dgemm_gcc

dgemm_icc: dgemm.c
	icc $(ICC_OPT) dgemm.c $(MKL_FLAGS) -lm -o dgemm_icc
