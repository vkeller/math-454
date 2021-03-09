all: a.out

a.out: main.cpp image-mat.o image-mat.h
	g++ main.cpp image-mat.o \
	  -I/usr/local/cuda/include \
	  -L/usr/local/cuda/lib64 \
	  -lcublas \
	  -lcudart

image-mat.o: image-mat.cu
	nvcc -c image-mat.cu  -arch sm_35
clean:
	rm *.o a.out
