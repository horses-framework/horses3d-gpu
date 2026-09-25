#################################################################################
##	Makefile for the WallShearStress unit test.
##
##	Requires the Navier-Stokes libraries to be built first, e.g. from ./Solver:
##	   make ns COMPILER=gfortran COMM=SERIAL ENABLE_THREADS=NO
##
##	Usage:
##	   make COMPILER=gfortran     (build)
##	   make run COMPILER=gfortran (build and run)
#################################################################################
HOME_DIR ?= ../../..
COMPILER ?= gfortran
PROG     := TestWallShearStress

ifeq ($(COMPILER),gfortran)
   FC     = gfortran
   FFLAGS = -cpp -ffree-line-length-0 -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
   MOD_FLAG = -J
else ifeq ($(COMPILER),ifort)
   FC     = ifort
   FFLAGS = -fpp -O0 -g -check all -traceback -fpe0
   MOD_FLAG = -module
else ifeq ($(COMPILER),nvfortran)
   FC     = nvfortran
   FFLAGS = -cpp -O0 -g
   MOD_FLAG = -module
endif

INCLUDE = -I$(HOME_DIR)/include -I$(HOME_DIR)/include_ns
LIBS    = -L$(HOME_DIR)/lib -lphysics_ns -lphysicsns -lio -lspectral -lmpiutils -lfoundation -lftobject

.DEFAULT_GOAL := all
all: ./bin/$(PROG)

./bin/$(PROG): ./src/$(PROG).f90
	@mkdir -p ./bin ./build
	$(FC) $(FFLAGS) $(INCLUDE) $(MOD_FLAG) ./build -o $@ $< $(LIBS)

run: ./bin/$(PROG)
	./bin/$(PROG)

clean:
	@rm -rf ./bin ./build

.PHONY: all run clean
