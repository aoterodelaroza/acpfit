## ifort with mkl
#FC = ifort
#FCFLAGS = -g -CU -C -traceback -fpe0 -debug -openmp
#LDFLAGS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -openmp
#FCFLAGS = -g -CU -C -traceback -fpe0 -debug
#LDFLAGS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
FCFLAGS = -O3 -openmp
LDFLAGS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -openmp

#### gfortran with static blas and lapack
FC = gfortran
#FCFLAGS = -O3 -ffree-line-length-none 
FCFLAGS = -g -fbounds-check -Wall -Wunused-parameter -ffpe-trap=invalid -fbacktrace -fdump-core -ffree-line-length-none
#LIBS=liblapack.a libblas.a 
LIBS=-llapack -lblas
LDFLAGS = 
#LDFLAGS =

BINS=acpfit
OBJS=files.o global.o types.o tools_io.o calc.o
INCLUDE=

%.o: %.f90
	$(FC) -c $(FCFLAGS) $(INCLUDE) -o $@ $<

%.o: %.f
	$(FC) -c $(FCFLAGS) $(INCLUDE) -o $@ $<

%.mod: %.o
	@if [ ! -f $@ ]; then rm $< ; $(MAKE) $< ; fi

# General targets

all: $(BINS)

clean:
	rm -f core *.mod *.o 

mrproper:
	rm -f core *.mod *.o $(BINS)

acpfit: $(OBJS) acpfit.o
	$(FC) -o acpfit $(LDFLAGS) $(OBJS) acpfit.o $(LIBS)

# Object dependencies
acpfit.o : calc.mod
acpfit.o : files.mod
files.o calc.o acpfit.o : global.mod
types.o tools.o global.o files.o calc.o acpfit.o : tools_io.mod
global.o calc.o acpfit.o : types.mod
