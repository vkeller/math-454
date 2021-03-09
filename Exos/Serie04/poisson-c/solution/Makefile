CC=gcc
LD=${CC}
CFLAGS+=-Wall -Werror -pedantic -O3 -fopenmp
LDFLAGS+=-lm ${CFLAGS}

OBJS=poisson.o io_binary.o timing.o simulation.o

all: poisson

poisson: $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)

clean:
	rm -f poisson *.o *~
