CXX=g++
LD=${CXX}
CXXFLAGS+=-Wall -Wextra -Werror -pedantic -std=c++11
LDFLAGS+=-lm $(CXXFLAGS)

OBJS=poisson.o simulation.o double_buffer.o grid.o dumpers.o

all: poisson

poisson: $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)

clean:
	rm -f hello poisson *.o *~
