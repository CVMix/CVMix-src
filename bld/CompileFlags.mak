ifeq ($(FC),ftn)
  ifeq ($(PE_ENV),GNU)
    FC_TMP = gfortran
  endif
  ifeq ($(PE_ENV),PGI)
    FC_TMP = pgf90
  endif
  ifeq ($(PE_ENV),INTEL)
    FC_TMP = ifx
  endif
  ifeq ($(PE_ENV),PATHSCALE)
    FC_TMP = pathf95
  endif
  ifeq ($(PE_ENV),CRAY)
    UCASE = TRUE
    FCFLAGS = -O2 -f free -e m -J $(OBJ_DIR)
  endif
else
  FC_TMP = $(FC)
endif

ifeq ($(FC_TMP),gfortran)
  FCFLAGS = -g -Og -Wall -Wno-maybe-uninitialized -Werror -fbacktrace -ffpe-trap=zero,overflow -fcheck=bounds -J $(OBJ_DIR)
endif

ifeq ($(FC_TMP),pgf90)
  FCFLAGS = -O2 -Mfree -module $(OBJ_DIR)
endif

ifeq ($(FC_TMP),ifx)
  FCFLAGS = -O2 -free -module $(OBJ_DIR) -cpp -warn all -diag-error warn -nogen-interface -fp-model source
endif

ifeq ($(FC_TMP),xlf90)
  FCFLAGS = -O2 -I $(OBJ_DIR) -qmoddir=$(OBJ_DIR)
endif

ifeq ($(FC_TMP),nagfor)
  FCFLAGS = -O2 -free -I $(OBJ_DIR) -mdir $(OBJ_DIR) -kind=byte
endif

ifeq ($(FC_TMP),pathf95)
  UCASE = TRUE
  FCFLAGS = -O2 -freeform -module $(OBJ_DIR)
endif
