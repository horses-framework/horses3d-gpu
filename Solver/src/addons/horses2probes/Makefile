#####################################################################
## HORSES - a High-Order Spectral Element Solver
##
##	horses2probes Makefile
##
#####################################################################

########################################
# 1.- User defined parameters
########################################
INSTALL_DIR?=./
FFLAGS?=-cpp
AR?=ar
ARFLAGS?=cr
RANLIB?=ranlib
MOD_FLAG?=-J
PLATFORM?=LINUX
COMPILER?=gfortran
MODE?=RELEASE
COMM?=PARALLEL
MAKE=make
AR=ar
RANLIB=ranlib
PROG=./bin/horses2probes

# Directory of shared source files from horses2tecplot
H2T_DIR=../horses2tecplot

# Files compiled from the horses2tecplot directory (shared code)
H2T_OBJS= Storage \
	OutputVariables \
	sharedSpectralBasis \

# Files compiled from the horses2probes directory (new code)
LOCAL_OBJS= Horses2probes \
	main \

H2T_OBJS_BUILD=$(foreach obj,$(H2T_OBJS),./build/$(obj).o)
LOCAL_OBJS_BUILD=$(foreach obj,$(LOCAL_OBJS),./build/$(obj).o)
ALL_OBJS_BUILD=$(H2T_OBJS_BUILD) $(LOCAL_OBJS_BUILD)


##########################################
##	Include dependencies 		##
##########################################
-include ./make.deps


##########################################################
##		COMPILATION RULES 			##
##########################################################

.DEFAULT_GOAL := $(PROG)

# Explicit rules for local files (explicit rules always beat pattern rules)
./build/Horses2probes.o: ./Horses2probes.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) $(MACROS) $(EXTLIB_INC) -I./include $(INCLUDE) $(INCLUDE)_ns $(MACROS) -c $< -o $@ $(MOD_FLAG) ./include

./build/main.o: ./main.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) $(MACROS) $(EXTLIB_INC) -I./include $(INCLUDE) $(INCLUDE)_ns $(MACROS) -c $< -o $@ $(MOD_FLAG) ./include

# Pattern rule for shared files from horses2tecplot
./build/%.o: $(H2T_DIR)/%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) $(MACROS) $(EXTLIB_INC) -I./include $(INCLUDE) $(INCLUDE)_ns $(MACROS) -c $< -o $@ $(MOD_FLAG) ./include

$(PROG): header mkdirs $(ALL_OBJS_BUILD)
	@echo
	@echo
	@echo "---------------------------------"
	@echo "| Linking " $@ "|"
	@echo "---------------------------------"
	$(FC) $(FFLAGS) $(EXTLIB_INC) -o $(PROG) $(ALL_OBJS_BUILD) $(HORSES_LIBS_EXEC) $(LIBS)


install: $(PROG)
	cp -p $(PROG) $(INSTALL_DIR)/bin
	cp -p ./include/*.mod $(INSTALL_DIR)/include
	@echo
	@echo

##########################################
##	Extra auxiliary rules 		##
##########################################
clean: FORCE
	rm -f ./include/*.mod ./lib/*.a ./build/*.o
	rm -f $(PROG)

mkdirs: FORCE
	@mkdir -p ./build
	@mkdir -p ./include
	@mkdir -p ./bin
	@mkdir -p ./lib

header: FORCE
	@echo
	@echo "================================"
	@echo ">>   Building horses2probes"
	@echo "================================"
	@echo 
	@echo

## Dummy target to make sure a rule is executed
FORCE: