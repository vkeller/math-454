CXX=mpic++
LD=${CXX}
CXXFLAGS+=-Wall -Wextra -Werror -pedantic -std=c++11 -O3
LDFLAGS+=-lm $(CXXFLAGS)

OBJS=poisson.o grid.o double_buffer.o simulation.o  dumpers.o

all: poisson

poisson: $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)

clean:
	rm -f hello poisson *.o *~
