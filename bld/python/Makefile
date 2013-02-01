FC=gfortran
RM=rm -f
FFLAGS=-ffree-form -ffree-line-length-none
FINCLUDES=

all: python

debug: python-debug

objects:
	(cd ../../src/shared; make FC="$(FC)" RM="$(RM)" FFLAGS="$(FFLAGS)" FINCLUDES="$(FINCLUDES)" )

objects-debug:
	(cd ../../src/shared; make FC="$(FC)" RM="$(RM)" FFLAGS="$(FFLAGS) -g" FINCLUDES="$(FINCLUDES)" )

python: objects
	rm -f cvmix.pyf cvmix.so
	f2py -m cvmix -h cvmix.pyf cvmix_BL_driver-pointers-python.F90
	cat cvmix.pyf | sed "s/inout/inplace/g" > temp
	mv temp cvmix.pyf
	f2py -c cvmix.pyf cvmix_BL_driver-pointers-python.F90 ../../src/shared/*.F90 -I../../src/shared --f90flags="$(FFLAGS)" --f77flags="$(FFLAGS)"

python-debug: objects-debug
	rm -f cvmix.pyf cvmix.so
	f2py -m cvmix -h cvmix.pyf cvmix_BL_driver-pointers-python.F90
	cat cvmix.pyf | sed "s/inout/inplace/g" > temp
	mv temp cvmix.pyf
	f2py -c cvmix.pyf cvmix_BL_driver-pointers-python.F90 ../../src/shared/*.F90 -I../../src/shared --f90flags="$(FFLAGS) -g" --f77flags="$(FFLAGS) -g" --debug-capi --debug


clean:
	(cd ../../src/shared; make clean)
	rm -f cvmix.pyf cvmix.so

