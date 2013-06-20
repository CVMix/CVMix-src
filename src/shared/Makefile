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
MAKE_DEP = $(SRC_DIR)/makedep.py
DEP_FILE = $(OBJ_DIR)/shared_deps.d

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

ifeq ($(USE_DEPS),TRUE)
  include $(DEP_FILE)
endif

### TARGETS ###

.PHONY: clean recurse check depends

all: depends recurse

recurse: check
ifeq ($(NEED_ENV_FILE),FALSE)
	$(MAKE) -e -f $(SRC_DIR)/Makefile $(LIB_DIR)/libcvmix.a $(INCS) USE_DEPS=TRUE
endif

# Create all object and module files
# Note that .mod files need to be copied to INC_DIR if OBJ_DIR != INC_DIR

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.F90
	$(FC) $(FCFLAGS) -c $< -o $@

$(INC_DIR)/%.mod: $(OBJ_DIR)/%.mod
ifneq ($(INC_DIR),$(OBJ_DIR))
	cp $< $@
endif

### Combine into library
$(LIB_DIR)/libcvmix.a: $(KINDS_OBJ) $(OBJS)
	ar -ru $(LIB_DIR)/libcvmix.a $(OBJS)

# Utilities:
# 1) make dependency file
# 2) check ensures that FC has been specified
# 3) clean deletes object and mod files
depends: $(DEP_FILE)

$(DEP_FILE): $(MAKE_DEP) $(MODULES)
	$(MAKE_DEP) $(SRC_DIR) $(OBJ_DIR) $(DEP_FILE) $(UCASE)
	@echo "Generated dependencies!"

check:
	@$(if $(NOFC), echo "ERROR: you must specify FC (and it is recommended that \
     you specify FCFLAGS"; echo "as well)."; exit 1)
ifeq ($(NEED_ENV_FILE),TRUE)
	$(BLD_DIR)/cvmix_setup $(ENV_FILE); $(MAKE) -e -f $(SRC_DIR)/Makefile
endif

clean:
	/bin/rm -f $(LIB_DIR)/libcvmix.a $(OBJS) $(KINDS_OBJ) $(MODS) $(KINDS_MOD) $(DEP_FILE)

