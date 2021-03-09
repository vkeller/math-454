################################################################################
# Automatically-generated file. Do not edit!
################################################################################

ARCH		= intel
# bluegene

ifeq ($(ARCH),bluegene)
CXX		= mpixlcxx_r
LD		= $(CXX)

CXXFLAGS	= -DBLUEGENE -DWRITEINFILE
CXXOPTS		= -O3 -qarch=qp -qtune=qp
LDFLAGS		=  -qarch=qp -qtune=qp

INCLUDES	= -I"../include" -I/bgsys/ibm_essl/prod/opt/ibmmath/essl/5.1/include
LAPACKDIR	= /bgsys/local/lib64
LIBS		= -L$(LAPACKDIR) -lesslbg -L/opt/ibmcmp/xlf/bg/14.1/lib64 -lxlf90_r -lxlsmp -L/opt/ibmcmp/xlmass/bg/7.3/lib64 -lmass
endif


ifeq ($(ARCH),intel)
CXX		= mpiicpc
LD		= $(CXX)

CXXFLAGS	= -DWRITEINFILE
CXXOPTS		= -O3
LDFLAGS		= 

INCLUDES	= -I"../include"
LAPACKDIR	= /usr/scratch/intel/mkl/lib/intel64
LIBS		= -L$(LAPACKDIR) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
endif


ifeq ($(ARCH),intel-phi)
CXX             = mpiicpc
LD              = $(CXX)

CXXFLAGS        = -DWRITEINFILE
CXXOPTS         = -O3 -xMIC-AVX512
LDFLAGS         = -xMIC-AVX512

INCLUDES        = -I"../include"
LAPACKDIR       = /usr/scratch/intel/mkl/lib/intel64
LIBS            = -L$(LAPACKDIR) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
endif


EXEC		= MXM-Cxx

RM 		= rm

.SUFFIXES:	.cxx .cpp .o

.cpp.o:	
		@echo 'Building target: $@'
		@echo 'Invoking: C++ Compiler on architecture $(ARCH)'
		@echo 'with this line : $(CXX) $(CXXOPTS) $(CXXFLAGS) $(INCLUDES) '
		$(CXX) $(CXXOPTS) $(CXXFLAGS) $(INCLUDES) -c $<
		@echo 'Finished compiling target: $@'
		@echo ' '

.cxx.o: 
		@echo 'Building target: $@'
		@echo 'Invoking: XL C++ Compiler'
		$(CXX) $(CXXOPTS) $(CXXFLAGS) $(INCLUDES) -c $<
		@echo 'Finished compiling target: $@'
		@echo ' '



OBJS		= MXM-Cxx.o mxm.o mxv.o

# All Target
all: MXM

# Tool invocations
MXM:	$(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: XL C++ Executable Linker'
	$(LD) $(OBJS) $(LIBS) $(LDFLAGS) -o $(EXEC)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS) $(EXEC) *.lst *.o *.dat *~
	-@echo ' '
