.SUFFIXES: .F90 .o
extension = .F90

# Need CVMix root directory
# If using old version of make, pass in CVMIX_ROOT
ifeq ($(wildcard $(MAKEFILE_LIST)),)
  CVMIX_ROOT := ../..
else
  ThisMakefile=$(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))
  ifeq (/,$(findstring /,$(ThisMakefile)))
    ThisDir=$(shell x=$(ThisMakefile) && echo $${x%/*})/
  endif
  CVMIX_ROOT := $(realpath $(ThisDir)../..)
endif

SETUP_FILE = $(CVMIX_ROOT)/bld/cvmix_setup
ENV_FILE = $(CVMIX_ROOT)/bld/.CVMix_env
# Directories used in Makefile
# Check to this is standalone CVMix or not (look for $(SETUP_FILE))
ifneq ($(wildcard $(SETUP_FILE)),)
  SRC_DIR = $(CVMIX_ROOT)/src/shared
  BLD_DIR = $(CVMIX_ROOT)/bld
  FLAG_FILE = $(BLD_DIR)/CompileFlags.mak
  OBJ_DIR = $(BLD_DIR)/obj
  LIB_DIR = $(CVMIX_ROOT)/lib
  INC_DIR = $(CVMIX_ROOT)/include
  ifeq ($(wildcard $(ENV_FILE)),)
    NEED_ENV_FILE = TRUE
  else
    # Read in Compiler information
    include $(ENV_FILE)
    NEED_ENV_FILE = FALSE
  endif
  # Read in FCFLAGS
  include $(FLAG_FILE)
else
  SRC_DIR = .
  OBJ_DIR = .
  LIB_DIR = .
  INC_DIR = .
  FC = NONE
  NEED_ENV_FILE = FALSE
endif

ifeq ($(FC),NONE)
  NOFC = TRUE
endif

MODULES = cvmix_background.F90 \
          cvmix_convection.F90 \
          cvmix_ddiff.F90 \
          cvmix_kpp.F90 \
          cvmix_put_get.F90 \
          cvmix_shear.F90 \
          cvmix_tidal.F90

# Some compilers produce ALL_UPPER_CASE.mod files
ifeq ($(UCASE),TRUE)
  MODS_TMP = CVMIX_KINDS_AND_TYPES.F90 \
             $(shell echo $(MODULES) | tr '[a-z]' '[A-Z]')
else
  MODS_TMP = cvmix_kinds_and_types.F90 \
             $(MODULES)
endif
ifneq ($(OBJ_DIR),$(INC_DIR))
  INCS = $(addprefix $(INC_DIR)/,${MODS_TMP:.F90=.mod})
endif
MODS = $(addprefix $(OBJ_DIR)/,${MODS_TMP:.F90=.mod}) \
       $(INCS)
OBJS = $(addprefix $(OBJ_DIR)/,${MODULES:.F90=.o})
KINDS_OBJ = $(OBJ_DIR)/cvmix_kinds_and_types.o
ifeq ($(FC),$(filter $(FC),pgf90 ftn))
  OBJS += $(KINDS_OBJ)
endif

### TARGETS ###

.PHONY: clean recurse check

all: recurse

recurse: check
ifeq ($(NEED_ENV_FILE),FALSE)
	$(MAKE) -e -f $(SRC_DIR)/Makefile $(LIB_DIR)/libcvmix.a $(INCS)
endif

# Create all object and module files
# Note that .mod files need to be copied to INC_DIR

### cvmix_kinds_and_types.F90
$(OBJ_DIR)/cvmix_kinds_and_types.o: $(SRC_DIR)/cvmix_kinds_and_types.F90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/cvmix_kinds_and_types.F90 -o $@

$(INC_DIR)/cvmix_kinds_and_types.mod: $(SRC_DIR)/cvmix_kinds_and_types.F90 \
                                      $(OBJ_DIR)/cvmix_kinds_and_types.o
	cp $(OBJ_DIR)/cvmix_kinds_and_types.mod $@

# Upper-case for Cray and pathscale compilers
$(INC_DIR)/CVMIX_KINDS_AND_TYPES.mod: $(SRC_DIR)/cvmix_kinds_and_types.F90 \
                                      $(OBJ_DIR)/cvmix_kinds_and_types.o
	cp $(OBJ_DIR)/CVMIX_KINDS_AND_TYPES.mod $@

### cvmix_put_get.F90
$(OBJ_DIR)/cvmix_put_get.o: $(SRC_DIR)/cvmix_put_get.F90       \
                            $(OBJ_DIR)/cvmix_kinds_and_types.o
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/cvmix_put_get.F90 -o $@

$(INC_DIR)/cvmix_put_get.mod: $(SRC_DIR)/cvmix_put_get.F90       \
                              $(OBJ_DIR)/cvmix_put_get.o
	cp $(OBJ_DIR)/cvmix_put_get.mod $@

$(INC_DIR)/CVMIX_PUT_GET.mod: $(SRC_DIR)/cvmix_put_get.F90       \
                              $(OBJ_DIR)/cvmix_put_get.o
	cp $(OBJ_DIR)/CVMIX_PUT_GET.mod $@

### cvmix_background.F90
$(OBJ_DIR)/cvmix_background.o: $(SRC_DIR)/cvmix_background.F90    \
                               $(OBJ_DIR)/cvmix_kinds_and_types.o
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/cvmix_background.F90 -o $@

$(INC_DIR)/cvmix_background.mod: $(SRC_DIR)/cvmix_background.F90       \
                                 $(OBJ_DIR)/cvmix_background.o
	cp $(OBJ_DIR)/cvmix_background.mod $@

$(INC_DIR)/CVMIX_BACKGROUND.mod: $(SRC_DIR)/cvmix_background.F90       \
                                 $(OBJ_DIR)/cvmix_background.o
	cp $(OBJ_DIR)/CVMIX_BACKGROUND.mod $@

### cvmix_convection.F90
$(OBJ_DIR)/cvmix_convection.o: $(SRC_DIR)/cvmix_convection.F90    \
                               $(OBJ_DIR)/cvmix_kinds_and_types.o
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/cvmix_convection.F90 -o $@

$(INC_DIR)/cvmix_convection.mod: $(SRC_DIR)/cvmix_convection.F90       \
                                 $(OBJ_DIR)/cvmix_convection.o
	cp $(OBJ_DIR)/cvmix_convection.mod $@

$(INC_DIR)/CVMIX_CONVECTION.mod: $(SRC_DIR)/cvmix_convection.F90       \
                                 $(OBJ_DIR)/cvmix_convection.o
	cp $(OBJ_DIR)/CVMIX_CONVECTION.mod $@

### cvmix_shear.F90
$(OBJ_DIR)/cvmix_shear.o: $(SRC_DIR)/cvmix_shear.F90         \
                          $(OBJ_DIR)/cvmix_background.o      \
                          $(OBJ_DIR)/cvmix_kinds_and_types.o
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/cvmix_shear.F90 -o $@

$(INC_DIR)/cvmix_shear.mod: $(SRC_DIR)/cvmix_shear.F90       \
                            $(OBJ_DIR)/cvmix_shear.o
	cp $(OBJ_DIR)/cvmix_shear.mod $@

$(INC_DIR)/CVMIX_SHEAR.mod: $(SRC_DIR)/cvmix_shear.F90       \
                            $(OBJ_DIR)/cvmix_shear.o
	cp $(OBJ_DIR)/CVMIX_SHEAR.mod $@

### cvmix_tidal.F90
$(OBJ_DIR)/cvmix_tidal.o: $(SRC_DIR)/cvmix_tidal.F90         \
                          $(OBJ_DIR)/cvmix_kinds_and_types.o
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/cvmix_tidal.F90 -o $@

$(INC_DIR)/cvmix_tidal.mod: $(SRC_DIR)/cvmix_tidal.F90       \
                            $(OBJ_DIR)/cvmix_tidal.o
	cp $(OBJ_DIR)/cvmix_tidal.mod $@

$(INC_DIR)/CVMIX_TIDAL.mod: $(SRC_DIR)/cvmix_tidal.F90       \
                            $(OBJ_DIR)/cvmix_tidal.o
	cp $(OBJ_DIR)/CVMIX_TIDAL.mod $@

### cvmix_ddiff.F90
$(OBJ_DIR)/cvmix_ddiff.o: $(SRC_DIR)/cvmix_ddiff.F90         \
                          $(OBJ_DIR)/cvmix_kinds_and_types.o
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/cvmix_ddiff.F90 -o $@

$(INC_DIR)/cvmix_ddiff.mod: $(SRC_DIR)/cvmix_ddiff.F90       \
                            $(OBJ_DIR)/cvmix_ddiff.o
	cp $(OBJ_DIR)/cvmix_ddiff.mod $@

$(INC_DIR)/CVMIX_DDIFF.mod: $(SRC_DIR)/cvmix_ddiff.F90       \
                            $(OBJ_DIR)/cvmix_ddiff.o
	cp $(OBJ_DIR)/CVMIX_DDIFF.mod $@

### cvmix_kpp.F90
$(OBJ_DIR)/cvmix_kpp.o: $(SRC_DIR)/cvmix_kpp.F90           \
                        $(OBJ_DIR)/cvmix_kinds_and_types.o
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/cvmix_kpp.F90 -o $@

$(INC_DIR)/cvmix_kpp.mod: $(SRC_DIR)/cvmix_kpp.F90         \
                          $(OBJ_DIR)/cvmix_kpp.o
	cp $(OBJ_DIR)/cvmix_kpp.mod $@

$(INC_DIR)/CVMIX_KPP.mod: $(SRC_DIR)/cvmix_kpp.F90         \
                          $(OBJ_DIR)/cvmix_kpp.o
	cp $(OBJ_DIR)/CVMIX_KPP.mod $@

### Combine into library
$(LIB_DIR)/libcvmix.a: $(KINDS_OBJ) $(OBJS)
	ar -ru $(LIB_DIR)/libcvmix.a $(OBJS)

# Utilities:
# 1) check ensures that FC has been specified
# 2) clean deletes object and mod files
check:
	@$(if $(NOFC), echo "ERROR: you must specify FC (and it is recommended that \
     you specify FCFLAGS"; echo "as well)."; exit 1)
ifeq ($(NEED_ENV_FILE),TRUE)
	$(BLD_DIR)/cvmix_setup $(ENV_FILE); $(MAKE) -e -f $(SRC_DIR)/Makefile
endif

clean:
	/bin/rm -f $(LIB_DIR)/libcvmix.a $(OBJS) $(KINDS_OBJ) $(MODS) $(KINDS_MOD)

