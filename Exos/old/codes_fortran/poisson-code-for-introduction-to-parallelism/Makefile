CC=mpicc
LD=${CC}
CFLAGS+=-Wall -Werror -g
LDFLAGS+=-lm

OBJS=poisson.o io_binary_mpi_io.o

all: poisson

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

poisson: $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)

hello : hello.o
	$(CC) $(CFLAGS) -o $@ $? $(LDFLAGS)

clean:
	rm -f hello poisson *.o *~
