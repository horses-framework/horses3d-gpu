#################################################################################
##	Makefile for the WallShearStress unit test.
##
##	It includes the solver Makefile, so it uses the same compiler, flags and
##	external libraries (COMPILER, MODE, COMM, WITH_METIS, EXTRA_FFLAGS, ...)
##	that were used to build HORSES3D. Only the equation-independent libraries
##	and the Navier-Stokes physics library (libphysicsns, built by any of
##	make ns / ins / mu) are needed, so it works after any solver build.
##
##	Usage (from this folder, after building the solver):
##	   make unittest     COMPILER=gfortran COMM=SEQUENTIAL   (build)
##	   make run-unittest COMPILER=gfortran COMM=SEQUENTIAL   (build and run)
#################################################################################
HOME_DIR ?= ../../..
include $(HOME_DIR)/Makefile

PROG_UNITTEST := TestWallShearStress
LIBS_UNITTEST := -L$(HOME_DIR)/lib -lphysicsns -lio -lspectral -lmpiutils -lfoundation -lftobject -lml
# Runtime search path for shared external libraries (e.g. HDF5 when WITH_HDF5=YES)
RPATH_UNITTEST := -Wl,-rpath,$(abspath $(HOME_DIR)/lib)
ifdef HDF5_DIR
   RPATH_UNITTEST += -Wl,-rpath,$(HDF5_DIR)/lib
endif

.DEFAULT_GOAL := unittest
unittest: ./bin/$(PROG_UNITTEST)

./bin/$(PROG_UNITTEST): ./src/$(PROG_UNITTEST).f90 FORCE
	@mkdir -p ./bin ./build
	$(FC) $(FFLAGS) $(MACROS) $(INCLUDE) $(EXTLIB_INC) $(MOD_FLAG) ./build -o $@ $< $(LIBS_UNITTEST) $(LIBS) $(RPATH_UNITTEST)

run-unittest: unittest
	./bin/$(PROG_UNITTEST)

.PHONY: unittest run-unittest
