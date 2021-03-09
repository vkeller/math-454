#OPTIM+=-O3 -march=native
DEBUG+=-g -O1
CXX=g++
CC=g++
LD=${CXX}
CXXFLAGS+=$(OPTIM) $(DEBUG) -Wall -Wextra -std=c++11
LDFLAGS+=$(OPTIM) $(DEBUG) $(CXXFLAGS) -lm

all: read write

pi: read.o write.o

clean:
	rm -f read write *.o *~
