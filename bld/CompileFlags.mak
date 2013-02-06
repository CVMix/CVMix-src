ifeq ($(FC),gfortran)
  COMPILE_FLAGS = -Wall -Werror -cpp -O2 -J $(OBJ_DIR)
endif

ifeq ($(FC),pgf90)
  COMPILE_FLAGS = -O2 -module $(OBJ_DIR)
endif

ifeq ($(FC),ifort)
  COMPILE_FLAGS = -diag-error warn,remark -warn all -cpp -O2 -module $(OBJ_DIR)
endif

ifeq ($(FC),xlf90)
  COMPILE_FLAGS = -O2 -I $(OBJ_DIR) -qmoddir=$(OBJ_DIR)
endif

ifeq ($(FC),nagfor)
  COMPILE_FLAGS = -kind=byte -free -O2 -I $(OBJ_DIR) -mdir $(OBJ_DIR)
endif
