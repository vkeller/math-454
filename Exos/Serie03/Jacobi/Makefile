CC=gcc
FLAGS  = -std=gnu99 -g -O3 -march=native -lm 
#
#FLAGS += -fopenmp 
#FLAGS += -DPAPI -I/ssoft/papi/5.4.0/gnu/include/ /ssoft/papi/5.4.0/gnu/lib/libpapi.a
#FLAGS +=  -I/apps/dora/iaca/iaca-lin32/include/
#
all: jacobi-naive jacobi-sse jacobi-avx jacobi-avx-peel jacobi-avx-block
#
jacobi-naive: jacobi.o jacobi-naive.o
	gcc $(FLAGS) $^ -o $@ 
#
jacobi-sse: jacobi.o jacobi-sse.o
	gcc $(FLAGS) $^ -o $@ 
#
jacobi-avx: jacobi.o jacobi-avx.o
	gcc $(FLAGS) $^ -o $@ 
#
jacobi-avx-peel: jacobi.o jacobi-avx-peel.o
	gcc $(FLAGS) $^ -o $@ 
#
jacobi-avx-block: jacobi.o jacobi-avx-block.o
	gcc $(FLAGS) $^ -o $@ 
#
%.o: %.c
	$(CC) -c $(FLAGS) $<
	$(CC) -S $(FLAGS) $< -o $@.s 
#
clean:
	rm -f *.s *.lst *.o jacobi-avx jacobi-avx-peel jacobi-avx-block jacobi-sse jacobi-naive
