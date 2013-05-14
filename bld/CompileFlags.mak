ifeq ($(FC),ftn)
  ifeq ($(PE_ENV),GNU)
    FC_TMP = gfortran
  endif
  ifeq ($(PE_ENV),PGI)
    FC_TMP = pgf90
  endif
  ifeq ($(PE_ENV),INTEL)
    FC_TMP = ifort
  endif
else
  FC_TMP = $(FC)
endif

ifeq ($(FC_TMP),gfortran)
  FCFLAGS = -Wall -Werror -cpp -O2 -J $(OBJ_DIR)
endif

ifeq ($(FC_TMP),pgf90)
  FCFLAGS = -O2 -module $(OBJ_DIR)
endif

ifeq ($(FC_TMP),ifort)
  FCFLAGS = -diag-error warn,remark -warn all -nogen-interface -cpp -O2 -module $(OBJ_DIR)
endif

ifeq ($(FC_TMP),xlf90)
  FCFLAGS = -O2 -I $(OBJ_DIR) -qmoddir=$(OBJ_DIR)
endif

ifeq ($(FC_TMP),nagfor)
  FCFLAGS = -kind=byte -free -O2 -I $(OBJ_DIR) -mdir $(OBJ_DIR)
endif
