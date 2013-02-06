.SUFFIXES: .F90 .o

# Need CVMix root directory
# If using old version of make, pass in CVMIX_ROOT
ifeq ($(wildcard $(MAKEFILE_LIST)),)
  CVMIX_ROOT = ../..
else
  ThisMakefile=$(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))
  ifeq (/,$(findstring /,$(ThisMakefile)))
    ThisDir=$(shell x=$(ThisMakefile) && echo $${x%/*})/
  endif
  CVMIX_ROOT = $(realpath $(ThisDir)../..)
endif

# Directory to build in
BLD = $(CVMIX_ROOT)/bld
OBJ = $(BLD)/obj
DIR = .

# Read in Compiler information
include $(BLD)/.CVMix_env
# Read in COMPILE_FLAGS
include $(BLD)/CompileFlags.mak

OBJS = $(OBJ)/cvmix_background.o \
	     $(OBJ)/cvmix_convection.o \
	     $(OBJ)/cvmix_ddiff.o \
	     $(OBJ)/cvmix_kinds_and_types.o \
	     $(OBJ)/cvmix_put_get.o \
	     $(OBJ)/cvmix_shear.o \
	     $(OBJ)/cvmix_tidal.o

all: libcvmix.a

# Create all object and module files
$(OBJ)/cvmix_kinds_and_types.o:
	$(FC) $(COMPILE_FLAGS) -c $(DIR)/cvmix_kinds_and_types.F90 -o $@

$(OBJ)/cvmix_put_get.o: $(OBJ)/cvmix_kinds_and_types.o
	$(FC) $(COMPILE_FLAGS) -c $(DIR)/cvmix_put_get.F90 -o $@

$(OBJ)/cvmix_background.o: $(OBJ)/cvmix_put_get.o         \
                          $(OBJ)/cvmix_kinds_and_types.o
	$(FC) $(COMPILE_FLAGS) -c $(DIR)/cvmix_background.F90 -o $@

$(OBJ)/cvmix_convection.o: $(OBJ)/cvmix_kinds_and_types.o
	$(FC) $(COMPILE_FLAGS) -c $(DIR)/cvmix_convection.F90 -o $@

$(OBJ)/cvmix_shear.o: $(OBJ)/cvmix_kinds_and_types.o
	$(FC) $(COMPILE_FLAGS) -c $(DIR)/cvmix_shear.F90 -o $@

$(OBJ)/cvmix_tidal.o: $(OBJ)/cvmix_kinds_and_types.o
	$(FC) $(COMPILE_FLAGS) -c $(DIR)/cvmix_tidal.F90 -o $@

$(OBJ)/cvmix_ddiff.o: $(OBJ)/cvmix_kinds_and_types.o
	$(FC) $(COMPILE_FLAGS) -c $(DIR)/cvmix_ddiff.F90 -o $@

libcvmix.a: $(OBJS)
	ar -ru $(DIR)/libcvmix.a $(OBJS)

clean:
	/bin/rm -f $(DIR)/*.a

