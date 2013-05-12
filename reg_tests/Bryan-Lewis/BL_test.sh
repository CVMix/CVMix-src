#!/bin/bash

##################
# Usage Subroutine
##################

usage() {
  echo './BL_test.sh [-mc|--memcopy] [-nc|--netcdf]'
  echo ''
  echo 'By default, this script allocates memory in the driver and the CVMix datatypes point back to memory. Output is ascii.'
  echo '  -h     | --help     Display this message'
  echo '  -mc    | --memcopy  Allocate memory in the CVMix datatypes rather than point'
  echo '  -nc    | --netcdf   Output to a netcdf file (data.nc) rather than three ascii files'
  echo '  -clean | --clean    Delete all output data (data.nc and any .out files)'
}

#############
# Main Script
#############

DRIVER=pointers
while [ $# -gt 0 ]; do
  case $1 in
    -h|--help)
      usage
      exit 0
      ;;
    -nc|--netcdf)
      USE_NETCDF=netcdf
      ;;
    -mc|--memcopy)
      DRIVER=mem_copy
      ;;
    -clean|--clean)
      echo "rm -f data.out col1.out col2.out data.nc"
      rm -f data.out col1.out col2.out data.nc
      exit
      ;;
    * )
      echo "ERROR: improper use of script. Run ./BL_test.sh -h for usage."
      exit 1
      ;;
  esac
  shift
done

CVMix=$PWD/../..
make -f $CVMix/src/Makefile ${USE_NETCDF}
# Note: if make error, include CVMIX_ROOT as below
#make -f $CVMix/src/Makefile ${USE_NETCDF} CVMIX_ROOT=$CVMix
if [ $? != 0 ]; then
  echo "Build error!"
  exit 1
fi

if [ "$DRIVER" == "mem_copy" ]; then
  $CVMix/bin/cvmix < input_memcopy.nl
else
  $CVMix/bin/cvmix < input_pointer.nl
fi

if [ $? != 0 ]; then
  echo "Error in execution!"
  exit 1
fi

if [ "${USE_NETCDF}" == "netcdf" ]; then
  if [ -f data.nc ]; then
    ncdump -v diff data.nc
  else
    echo "Execution error!"
    exit 1
  fi
else
  if [ -f data.out ]; then
    cat data.out
  else
    echo "Execution error!"
    exit 1
  fi
fi
