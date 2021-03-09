CC=gcc
LD=${CC}
DEBUG=-g
CFLAGS+=-Wall -Wextra -pedantic -O3  $(DEBUG)
LDFLAGS+=-lm -lm $(DEBUG)

OBJS_NBODY=nbody.o nbody_bruteforce.o reader.o nbody_barneshut.o

all: nbody-code

nbody-code: $(OBJS_NBODY)
	$(LD) $(OBJS_NBODY) $(LDFLAGS) -o $@

clean:
	rm -Rf nbody-code *.o *~  gmon.out
