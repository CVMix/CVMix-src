.SUFFIXES: .F90 .o

OBJS = cvmix_background.o \
	   cvmix_convection.o \
	   cvmix_ddiff.o \
	   cvmix_kinds_and_types.o \
	   cvmix_put_get.o \
	   cvmix_shear.o \
	   cvmix_tidal.o \

all: $(OBJS)

cvmix_kinds_and_types.o:

cvmix_put_get.o: cvmix_kinds_and_types.o

cvmix_background.o: cvmix_put_get.o

cvmix_convection.o: cvmix_kinds_and_types.o

cvmix_ddiff.o: cvmix_kinds_and_types.o

cvmix_shear.o: cvmix_kinds_and_types.o

cvmix_tidal.o: cvmix_kinds_and_types.o

clean:
	$(RM) *.o *.mod

.F90.o:
	$(FC) $(FFLAGS) -c $< $(FCINCLUDES)

