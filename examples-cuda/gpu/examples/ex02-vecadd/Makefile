# Make file for vectorAdd example

# Location of the CUDA Toolkit binaries and libraries
CUDA_PATH       ?= /usr/local/apps/cuda-5.0
CUDA_INC_PATH   ?= $(CUDA_PATH)/include
CUDA_BIN_PATH   ?= $(CUDA_PATH)/bin
CUDA_LIB_PATH  ?= $(CUDA_PATH)/lib64

# Common binaries
NVCC            ?= $(CUDA_BIN_PATH)/nvcc
GCC             ?= g++

# CUDA code generation flags
GENCODE_FLAGS    := -gencode arch=compute_20,code=sm_20

# OS-specific build flags
LDFLAGS   := -L$(CUDA_LIB_PATH) -lcudart
CCFLAGS   := -m64

# OS-architecture specific flags
NVCCFLAGS := -m64

# Target rules
all: build

build: vectorAdd

vectorAdd.o: vectorAdd.cu
	$(NVCC) $(NVCCFLAGS) -I$(CUDA_INC_PATH) $(GENCODE_FLAGS) -o $@ -c $<

vectorAdd: vectorAdd.o
	$(GCC) $(CCFLAGS) -o $@ $+ $(LDFLAGS)

run: build
	./vectorAdd

clean:
	rm -f vectorAdd vectorAdd.o 
