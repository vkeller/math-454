CC=mpicc
CFLAGS=-Wall -Werror -g
LDFLAGS=-lm

EXECUTABLES=hello ping ping_correct iping exchange exchange_send_recv exchange_send_recv_replace

all: $(EXECUTABLES)

hello: hello.c
	$(CC) $< $(CFLAGS) -o $@

ping: ping.c
	$(CC) $< $(CFLAGS) -o $@

ping_correct: ping_correct.c
	$(CC) $< $(CFLAGS) -o $@

iping: iping.c
	$(CC) $< $(CFLAGS) -o $@

exchange: exchange.c
	$(CC) $< $(CFLAGS) -o $@

exchange_send_recv: exchange_send_recv.c
	$(CC) $< $(CFLAGS) -o $@

exchange_send_recv_replace: exchange_send_recv_replace.c
	$(CC) $< $(CFLAGS) -o $@

clean:
	rm -f $(EXECUTABLES) timing.o
