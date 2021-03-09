OPTIM+=-O3

CXX=mpicxx
CC=mpicxx
LD=${CXX}
CXXFLAGS+=-Wall -Wextra -std=c++11 $(OPTIM)
LDFLAGS+=$(OPTIM) -lm

EXECS=pi_io_at pi_p2p_permanent_ring pi_p2p_derived_type

all: clean $(EXECS)

clean:
	rm -f $(EXECS) *.o *~
