# Makefile to create the 3D cosmic shear code

#compilers 
CC=/OPT/INTEL/bin/icc
FC=/OPT/INTEL/bin/ifort

#for cc (add -pg to profile using gprof)
#CFLAGS= -m64 -lm -lc -O3 -m64 -fPIC -O3 -ftree-vectorize -funroll-loops -w -fopenmp
#FFLAGS= -m64 -lm -O3 -m64 -fPIC -O3 -ftree-vectorize -funroll-loops -w -fopenmp

#for icc 
CFLAGS= -fast -openmp -L/OPT/INTEL/lib/intel64/ -lifcore -Wall
FFLAGS= -fast -openmp -nologo -cxxlib -nofor-main -L/OPT/INTEL/lib/intel64/ -lifcore -Wall

#platform
PLAT = _LINUX

# posix threads library name (not needed on some machines)
THREADLIB=-lpthread

#general library path
LIBDIR=${PWD}/../..

#LAPACK library
LAPATH=-L${LIBDIR}/libs/CLAPACK
LAINCLUDE=-I${LIBDIR}/libs/CLAPACK

#GNU science library 
GSLPATH=-L${LIBDIR}/libs/gsl-1.12/lib -lgsl -lgslcblas
GSLINCLUDE=-I${LIBDIR}/libs/gsl-1.12

# general include path
INCLUDE=-I/OPT/INTEL/include/
TDINCLUDE=-I${PWD}/../include/

# general library path
LIBRARY=-L/OPT/INTEL/lib/intel64/

# cfitsio library path, library name 
FITSPATH=-L${LIBDIR}/libs/cfitsio
FITSLIB=-lcfitsio
FITSINCLUDE=-I${LIBDIR}/libs/cfitsio

#fftw 
FFTWINCLUDE=-I${LIBDIR}/libs/fftw-3.2.1/include
FFTWLIB=-L${LIBDIR}/libs/fftw-3.2.1/lib/ -lfftw3

matrixobjects=\
write3Dmatrix_fits.o \
write2Dmatrix_fits.o \

besselobjects=\
smk_shbessel.o

cosmoobjects= \
commonfunctions.o 

# targets
smk_transform_data: smk_transform_data_v1.o  $(matrixobjects) $(besselobjects) $(cosmoobjects)
	$(CC) -o $@  smk_transform_data_v1.o \
	$(matrixobjects) $(besselobjects) $(cosmoobjects) \
	$(FORTRANLINK) $(THREADLIB) \
	$(LAPATH) -ltmg -llapack -lblas $(LAPATH)/F2CLIBS -lI77 -lF77 \
	$(GSLPATH) $(GSLINCLUDE) \
	$(CFLAGS) $(INCLUDE) $(TDINCLUDE) $(LIBRARY) \
	$(FITSPATH) $(FITSLIB) \
	$(FFTWINCLUDE) $(FFTWLIB) \
	$(OPTIM) $(WARNINGS) $(ARCHFLAGS) $(FPIC)

all: smk_transform_data

# general compilation rules 
.SUFFIXES:  .c   .o
.c.o:
	$(CC) -c $< -O3 -o $@ \
	$(CFLAGS) \
	$(FORTRANLINK) $(THREADLIB) \
	$(LAPATH) -ltmg -llapack -lblas $(LAPATH)/F2CLIBS -lI77 -lF77 \
	$(GSLPATH) $(GSLINCLUDE) \
	$(CFLAGS) $(INCLUDE) $(TDINCLUDE) $(LIBRARY) $(FITSINCLUDE) \
	$(FFTWINCLUDE) $(FFTWLIB) \
	$(OPTIM) $(WARNINGS) $(ARCHFLAGS) $(FPIC)

# general compilation for fortran 90 

.SUFFIXES:  .f90 .o
.f90.o:
	$(FC) -O3  $(FFLAGS) -c $< 

.SUFFIXES:  .for .o
.for.o:
	$(FC) -O3  $(FFLAGS) -c $< 

.SUFFIXES:  .f .o
.f.o:
	$(FC) -O3  $(FFLAGS) -c $<  $(FFLAGS)

tidy:	
	rm -f *.o 

clean:
	rm -f *.o *.mod smk_transform_data *~

