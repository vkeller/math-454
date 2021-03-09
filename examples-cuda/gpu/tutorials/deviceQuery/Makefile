# Makefile for Device Query program

# Location of the CUDA Toolkit binaries and libraries
CUDA_PATH       ?= $(SCC_CUDA_DIR)/install
CUDA_INC_PATH   ?= $(CUDA_PATH)/include
CUDA_BIN_PATH   ?= $(CUDA_PATH)/bin
CUDA_LIB_PATH  ?= $(CUDA_PATH)/lib64

# Common binaries
NVCC            ?= $(CUDA_BIN_PATH)/nvcc
GCC             ?= g++

# CUDA code generation flags
GENCODE_FLAGS    := -gencode arch=compute_20,code=sm_20

# OS-specific build flags
LDFLAGS   := -L$(CUDA_LIB_PATH) -lcuda -lcudart
CCFLAGS   := -m64

# OS-architecture specific flags
NVCCFLAGS := -m64

# Common includes and paths for CUDA
INCLUDES      := -I$(CUDA_INC_PATH) -I. 

# Target rules
all: build

build: deviceQuery

deviceQuery.o: deviceQuery.cpp
	$(GCC) $(INCLUDES) -o $@ -c $<

deviceQuery: deviceQuery.o
	$(GCC) -o $@ $+ $(LDFLAGS)

run: build
	./deviceQuery

clean:
	rm -f deviceQuery deviceQuery.o 
