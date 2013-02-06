.SUFFIXES: .F90 .o
extension = .F90

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
INC_DIR = $(CVMIX_ROOT)/include

# Read in Compiler information
include $(BLD_DIR)/.CVMix_env
# Read in FCFLAGS
include $(BLD_DIR)/CompileFlags.mak

SRC = cvmix_background.F90 \
	    cvmix_convection.F90 \
	    cvmix_ddiff.F90 \
	    cvmix_put_get.F90 \
	    cvmix_shear.F90 \
	    cvmix_tidal.F90

OBJS = $(addprefix $(OBJ_DIR)/,${SRC:.F90=.o})
INCS = $(addprefix $(INC_DIR)/,${SRC:.F90=.mod})

### TARGETS ###

.PHONY: clean recurse

all: recurse

recurse:
	$(MAKE) $(LIB_DIR)/libcvmix.a $(INCS) $(INC_DIR)/cvmix_kinds_and_types.mod

# Create all object and module files
# Note that .mod files need to be copied to INC_DIR

### cvmix_kinds_and_types.F90
$(OBJ_DIR)/cvmix_kinds_and_types.o: $(SRC_DIR)/cvmix_kinds_and_types.F90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/cvmix_kinds_and_types.F90 -o $@

$(INC_DIR)/cvmix_kinds_and_types.mod: $(SRC_DIR)/cvmix_kinds_and_types.F90 \
                          $(OBJ_DIR)/cvmix_kinds_and_types.o
	cp $(OBJ_DIR)/cvmix_kinds_and_types.mod $@

### cvmix_put_get.F90
$(OBJ_DIR)/cvmix_put_get.o: $(SRC_DIR)/cvmix_put_get.F90       \
                            $(OBJ_DIR)/cvmix_kinds_and_types.o
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/cvmix_put_get.F90 -o $@

$(INC_DIR)/cvmix_put_get.mod: $(SRC_DIR)/cvmix_put_get.F90       \
                              $(OBJ_DIR)/cvmix_put_get.o
	cp $(OBJ_DIR)/cvmix_put_get.mod $@

### cvmix_background.F90
$(OBJ_DIR)/cvmix_background.o: $(SRC_DIR)/cvmix_background.F90    \
                               $(OBJ_DIR)/cvmix_put_get.o         \
                               $(OBJ_DIR)/cvmix_kinds_and_types.o
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/cvmix_background.F90 -o $@

$(INC_DIR)/cvmix_background.mod: $(SRC_DIR)/cvmix_background.F90       \
                              $(OBJ_DIR)/cvmix_background.o
	cp $(OBJ_DIR)/cvmix_background.mod $@

### cvmix_convection.F90
$(OBJ_DIR)/cvmix_convection.o: $(SRC_DIR)/cvmix_convection.F90    \
                               $(OBJ_DIR)/cvmix_kinds_and_types.o
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/cvmix_convection.F90 -o $@

$(INC_DIR)/cvmix_convection.mod: $(SRC_DIR)/cvmix_convection.F90       \
                              $(OBJ_DIR)/cvmix_convection.o
	cp $(OBJ_DIR)/cvmix_convection.mod $@

### cvmix_shear.F90
$(OBJ_DIR)/cvmix_shear.o: $(SRC_DIR)/cvmix_shear.F90         \
                          $(OBJ_DIR)/cvmix_kinds_and_types.o
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/cvmix_shear.F90 -o $@

$(INC_DIR)/cvmix_shear.mod: $(SRC_DIR)/cvmix_shear.F90       \
                              $(OBJ_DIR)/cvmix_shear.o
	cp $(OBJ_DIR)/cvmix_shear.mod $@

### cvmix_tidal.F90
$(OBJ_DIR)/cvmix_tidal.o: $(SRC_DIR)/cvmix_tidal.F90         \
                          $(OBJ_DIR)/cvmix_kinds_and_types.o
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/cvmix_tidal.F90 -o $@

$(INC_DIR)/cvmix_tidal.mod: $(SRC_DIR)/cvmix_tidal.F90       \
                              $(OBJ_DIR)/cvmix_tidal.o
	cp $(OBJ_DIR)/cvmix_tidal.mod $@

### cvmix_ddiff.F90
$(OBJ_DIR)/cvmix_ddiff.o: $(SRC_DIR)/cvmix_ddiff.F90         \
                          $(OBJ_DIR)/cvmix_kinds_and_types.o
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/cvmix_ddiff.F90 -o $@

$(INC_DIR)/cvmix_ddiff.mod: $(SRC_DIR)/cvmix_ddiff.F90       \
                              $(OBJ_DIR)/cvmix_ddiff.o
	cp $(OBJ_DIR)/cvmix_ddiff.mod $@

### Combine into library
$(LIB_DIR)/libcvmix.a: $(OBJS)
	ar -ru $(LIB_DIR)/libcvmix.a $(OBJS)

clean:
	/bin/rm -f $(LIB_DIR)/libcvmix.a $(INC_DIR)/cvmix_kinds_and_types.mod

