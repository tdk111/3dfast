# Makefile to create the 3D cosmic shear code
# Thomas Kitching, 2015

#compilers 
CC=gcc
FC=gcc

#flags (add -pg to profile using gprof)
CFLAGS= -m64 -lm -lc -O3 
FFLAGS= -m64 -lm -O3 

#platform
PLAT = _LINUX

WARNINGS   = -Wall -Wuninitialized -pedantic
OPTIM      = -O3
FPIC	   = -fPIC
ARCHFLAGS  = -std=gnu9x

all:
	cd ${PWD}/src; make smk_transform_data
	cp ${PWD}/src/smk_transform_data ${PWD}/bin/

tidy:	
	rm -f *.o 

clean:
	cd ${PWD}/src; make clean
	rm -f ${PWD}/src/*.o ${PWD}/src/*.mod ${PWD}/bin/smk* *~ */*~

