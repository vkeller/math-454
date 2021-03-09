OPTIM+=-O3

CXX=mpicxx
CC=mpicxx
LD=${CXX}
CXXFLAGS+=-Wall -Wextra -std=c++11 $(OPTIM)
LDFLAGS+=$(OPTIM) -lm

EXECS=pi pi_p2p_ring

all: clean $(EXECS)

clean:
	rm -f $(EXECS) *.o *~
