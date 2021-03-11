CXX=g++
LD=${CXX}
CXXFLAGS+=
LDFLAGS+=
OBJS=hello.o greetings.o

hello: $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

clean:
	rm -f hello *.o *~
