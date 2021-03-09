CC = icc 

CFLAGS = -O3 -march=native 
CFLAGS += -Wall
CFLAGS += -qopenmp


all: stream 

stream: stream.c
	$(CC) $(CFLAGS) -qopt-streaming-stores=always stream.c -o stream

clean:
	rm -f stream *.o
