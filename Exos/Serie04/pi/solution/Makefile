OPTIM+=-O3 -march=native -fopenmp
CXX=g++
CC=g++
LD=${CXX}
CXXFLAGS+=-Wall -Wextra -std=c++11 $(OPTIM)
LDFLAGS+=-lm

EXE=pi_for_wrong pi_critical pi_critical_correct pi_reduction

all: clean $(EXE)

clean:
	rm -f $(EXE) *.o *~
