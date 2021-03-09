CC=gcc
CFLAGS=-Wall -Werror -g -fopenmp
LDFLAGS=-lm

EXECUTABLES=hello reduction scheduler

all: $(EXECUTABLES)

hello: hello.c
	$(CC) $< $(CFLAGS) -o $@

reduction: reduction.c
	$(CC) $< $(CFLAGS) -o $@

scheduler: scheduler.c timing.o
	$(CC) $< timing.o $(CFLAGS) -o $@ $(LDFLAGS)

clean:
	rm -f $(EXECUTABLES) timing.o
