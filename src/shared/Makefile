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

# Directories used in Makefile
SRC_DIR = $(CVMIX_ROOT)/src/shared
BLD_DIR = $(CVMIX_ROOT)/bld
OBJ_DIR = $(BLD_DIR)/obj
LIB_DIR = $(CVMIX_ROOT)/lib

# Read in Compiler information
include $(BLD_DIR)/.CVMix_env
# Read in COMPILE_FLAGS
include $(BLD_DIR)/CompileFlags.mak

OBJS = $(OBJ_DIR)/cvmix_background.o \
	     $(OBJ_DIR)/cvmix_convection.o \
	     $(OBJ_DIR)/cvmix_ddiff.o \
	     $(OBJ_DIR)/cvmix_put_get.o \
	     $(OBJ_DIR)/cvmix_shear.o \
	     $(OBJ_DIR)/cvmix_tidal.o

### TARGETS ###

.PHONY: clean recurse

all: recurse

recurse:
	$(MAKE) $(LIB_DIR)/libcvmix.a

# Create all object and module files
$(OBJ_DIR)/cvmix_kinds_and_types.o: $(SRC_DIR)/cvmix_kinds_and_types.F90
	$(FC) $(COMPILE_FLAGS) -c $(SRC_DIR)/cvmix_kinds_and_types.F90 -o $@

$(OBJ_DIR)/cvmix_put_get.o: $(SRC_DIR)/cvmix_put_get.F90       \
                            $(OBJ_DIR)/cvmix_kinds_and_types.o
	$(FC) $(COMPILE_FLAGS) -c $(SRC_DIR)/cvmix_put_get.F90 -o $@

$(OBJ_DIR)/cvmix_background.o: $(SRC_DIR)/cvmix_background.F90    \
                               $(OBJ_DIR)/cvmix_put_get.o         \
                               $(OBJ_DIR)/cvmix_kinds_and_types.o
	$(FC) $(COMPILE_FLAGS) -c $(SRC_DIR)/cvmix_background.F90 -o $@

$(OBJ_DIR)/cvmix_convection.o: $(SRC_DIR)/cvmix_convection.F90    \
                               $(OBJ_DIR)/cvmix_kinds_and_types.o
	$(FC) $(COMPILE_FLAGS) -c $(SRC_DIR)/cvmix_convection.F90 -o $@

$(OBJ_DIR)/cvmix_shear.o: $(SRC_DIR)/cvmix_shear.F90         \
                          $(OBJ_DIR)/cvmix_kinds_and_types.o
	$(FC) $(COMPILE_FLAGS) -c $(SRC_DIR)/cvmix_shear.F90 -o $@

$(OBJ_DIR)/cvmix_tidal.o: $(SRC_DIR)/cvmix_tidal.F90         \
                          $(OBJ_DIR)/cvmix_kinds_and_types.o
	$(FC) $(COMPILE_FLAGS) -c $(SRC_DIR)/cvmix_tidal.F90 -o $@

$(OBJ_DIR)/cvmix_ddiff.o: $(SRC_DIR)/cvmix_ddiff.F90         \
                          $(OBJ_DIR)/cvmix_kinds_and_types.o
	$(FC) $(COMPILE_FLAGS) -c $(SRC_DIR)/cvmix_ddiff.F90 -o $@

$(LIB_DIR)/libcvmix.a: $(OBJS)
	ar -ru $(LIB_DIR)/libcvmix.a $(OBJS)

clean:
	/bin/rm -f $(LIB_DIR)/*.a

