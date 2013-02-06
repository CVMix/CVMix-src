ifeq ($(FC),gfortran)
  COMPILE_FLAGS = -Wall -Werror -cpp -O2 -J $(OBJ)
endif

ifeq ($(FC),pgf90)
  COMPILE_FLAGS = -O2 -module $(OBJ)
endif

ifeq ($(FC),ifort)
  COMPILE_FLAGS = -diag-error warn,remark -warn all -cpp -O2 -module $(OBJ)
endif

ifeq ($(FC),xlf90)
  COMPILE_FLAGS = -O2 -I $(OBJ) -qmoddir=$(OBJ)
endif

ifeq ($(FC),nagfor)
  COMPILE_FLAGS = -kind=byte -free -O2 -I $(OBJ) -mdir $(OBJ)
endif
