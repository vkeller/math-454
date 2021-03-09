include ./make.inc


all: dgemm


HEADERS1 = utils.h 
SOURCES1 = dgemm.c  

OBJ1 = $(SOURCES1:.c=.o)

dgemm: $(OBJ1)
	$(CC) $(OBJ1) $(LDOPTS) $(OPTS) $(LIBDIR) $(LIB) $(INC) -o dgemm 

#.SUFFIXES: .c .cpp .cu .o

%.o : %.c
	echo $(OBJ)
	$(CC) $(OPTS) $(INC) -c $<


clean:
	echo $(OBJ1) dgemm 
	rm -f $(OBJ1) dgemm  


.PHONY: clean

