ifeq ($(FC),gfortran)
  FCFLAGS = -Wall -Werror -cpp -O2 -J $(OBJ_DIR)
endif

ifeq ($(FC),pgf90)
  FCFLAGS = -O2 -module $(OBJ_DIR)
endif

ifeq ($(FC),ifort)
  FCFLAGS = -diag-error warn,remark -warn all -nogen-interface -cpp -O2 -module $(OBJ_DIR)
endif

ifeq ($(FC),xlf90)
  FCFLAGS = -O2 -I $(OBJ_DIR) -qmoddir=$(OBJ_DIR)
endif

ifeq ($(FC),nagfor)
  FCFLAGS = -kind=byte -free -O2 -I $(OBJ_DIR) -mdir $(OBJ_DIR)
endif

ifeq ($(FC),ftn)
  FCFLAGS = -Mfree -O2 -module $(OBJ_DIR)
endif
