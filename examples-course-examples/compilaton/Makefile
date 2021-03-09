CC = gcc
CFLAGS = -O3 -ftree-vectorize
LDFLAGS = -lm

all: code

# This make file cut the compilation in 3 phases only for the purpose to show the phases

# compilation : assembly to machine code
code.o: code.s
	@echo ""; echo ""
	@echo " >> $(CC) -c $(CFLAGS) $< -o $@"
	@$(CC) -c $(CFLAGS) $< -o $@
	@objdump -d $@

# compilation : code to assembly code
code.s: code.i
	@echo ""; echo ""
	@echo " >> $(CC) -S $(CFLAGS) $< -o $@"
	@$(CC) -S $(CFLAGS) $< -o $@
	@cat $@

# precompilation
code.i: code.c
	@echo ""; echo ""
	@echo " >> $(CC) -E $(CFLAGS) $< -o $@"
	@$(CC) -E $(CFLAGS) $< -o $@
	@tail -n 12 $@

# linking
code: code.o
	@echo ""; echo ""
	@echo " >> $(CC) $(LDFLAGS) code.o -o $@"
	@$(CC) $(LDFLAGS) code.o -o $@

clean:
	rm -f *.o code code.i code.s
