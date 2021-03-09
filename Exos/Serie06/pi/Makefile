OPTIM+=-O3 -march=native
CXX=g++
CC=g++
LD=${CXX}
CXXFLAGS+=-Wall -Wextra -std=c++11 $(OPTIM)
LDFLAGS+=-lm
EXE=pi

all: clean $(EXE)

clean:
	rm -f $(EXE) *.o *~
