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
make -f $CVMix/bld/Makefile ${USE_NETCDF} VMIX_DRIVER=vmix_BL_driver-${DRIVER}.F90
# Note: if make error, include VMIX_ROOT as below
#make -f $CVMix/bld/Makefile ${USE_NETCDF} VMIX_DRIVER=vmix_BL_driver-${DRIVER}.F90 VMIX_ROOT=$CVMix
if [ $? != 0 ]; then
  echo "Build error!"
  exit 1
fi

$CVMix/bld/exe/cvmix < input.nl
if [ "${USE_NETCDF}" == "netcdf" ]; then
  ncdump -v diff data.nc
else
  cat data.out
fi
