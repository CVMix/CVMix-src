#!/bin/bash

##################
# Usage Subroutine
##################

usage() {
  echo './BL_test.sh [-mc|--memcopy] [-nc|--netcdf]'
  echo ''
  echo 'By default, this script allocates memory in the driver and the CVMix datatypes point back to memory. Output is ascii.'
  echo '  -h     | --help      Display this message'
  echo '  -mc    | --memcopy   Allocate memory in the CVMix datatypes rather than point'
  echo '  -nc    | --netcdf    Output to a netcdf file (data.nc) rather than three ascii files'
  echo '  -nb    | --no_build  Do not rebuild the executable'
  echo '  -clean | --clean     Delete all output data (data.nc and any .out files)'
}

#############
# Main Script
#############

DRIVER=pointers
CVMix=$( cd ../.. ; pwd )

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
    -nb|--nobuild)
      if [ ! -e $CVMix/bin/cvmix ]; then
        echo "ERROR: executable not found!"
        exit 1
      fi
      NO_BUILD=1
      ;;
    * )
      echo "ERROR: improper use of script. Run ./BL_test.sh -h for usage."
      exit 1
      ;;
  esac
  shift
done

if [ -z ${NO_BUILD} ]; then
  make -f $CVMix/src/Makefile CVMIX_ROOT=$CVMix ${USE_NETCDF}
  if [ $? != 0 ]; then
    echo "Build error!"
    exit 2
  fi
else
  if [ "`cat ../../bld/.netcdf_info`" == "YES" ]; then
    USE_NETCDF=netcdf
  fi  
fi

if [ "$DRIVER" == "mem_copy" ]; then
  $CVMix/bin/cvmix < input_memcopy.nl
else
  $CVMix/bin/cvmix < input_pointer.nl
fi

if [ $? != 0 ]; then
  echo "Error in execution!"
  exit 3
fi

if [ "${USE_NETCDF}" == "netcdf" ]; then
  if [ -f data.nc ]; then
    ncdump -v Tdiff data.nc
  else
    echo "ERROR: Output not found!"
    exit 4
  fi
else
  if [ -f data.out ]; then
    cat data.out
  else
    echo "ERROR: Output not found!"
    exit 4
  fi
fi
