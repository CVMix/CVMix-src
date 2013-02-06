#!/bin/bash

##################
# Usage Subroutine
##################

usage() {
  echo './double_diff-test.sh [-nc|--netcdf]'
  echo ''
  echo 'By default, output is ascii'
  echo '  -h     | --help     Display this message'
  echo '  -nc    | --netcdf   Output to a netcdf file (data.nc) rather than three ascii files'
  echo '  -clean | --clean    Delete all output data (data.nc and data.out files)'
}

#############
# Main Script
#############

while [ $# -gt 0 ]; do
  case $1 in
    -h|--help)
      usage
      exit 0
      ;;
    -nc|--netcdf)
      USE_NETCDF=netcdf
      ;;
    -clean|--clean)
      echo "rm -f data.out data.nc"
      rm -f data.out data.nc
      exit
      ;;
    * )
      echo "ERROR: improper use of script. Run ./double_diff_test.sh -h for usage."
      exit 1
      ;;
  esac
  shift
done

CVMix=$PWD/../..
make -f $CVMix/bld/Makefile ${USE_NETCDF} CVMIX_DRIVER=cvmix_ddiff_driver.F90
# Note: if make error, include CVMIX_ROOT as below
#make -f $CVMix/bld/Makefile ${USE_NETCDF} CVMIX_DRIVER=cvmix_ddiff_driver.F90 CVMIX_ROOT=$CVMix
if [ $? != 0 ]; then
  echo "Build error!"
  exit 1
fi

$CVMix/bin/cvmix < input.nl
if [ "${USE_NETCDF}" == "netcdf" ]; then
  ncdump -v visc data.nc
else
  cat data.out
fi
