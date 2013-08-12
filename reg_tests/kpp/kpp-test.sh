#!/bin/bash

##################
# Usage Subroutine
##################

usage() {
  echo './kpp-test.sh [-nc|--netcdf]'
  echo ''
  echo 'By default, output is ascii'
  echo '  -h     | --help     Display this message'
  echo '  -nc    | --netcdf   Output to a netcdf file (data.nc) rather than ascii file'
  echo '  -clean | --clean    Delete all output data (data.nc and data.out files)'
}

#############
# Main Script
#############

OUTPUT=TRUE
while [ $# -gt 0 ]; do
  case $1 in
    -h|--help)
      usage
      exit 0
      ;;
    -nc|--netcdf)
      USE_NETCDF=netcdf
      ;;
    -no_out|--no_output)
      OUTPUT=FALSE
      ;;
    -clean|--clean)
      cln_cmd="rm -f test*.out test*.nc"
      echo ${cln_cmd}
      ${cln_cmd}
      exit
      ;;
    * )
      echo "ERROR: improper use of script. Run ./kpp_test.sh -h for usage."
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

$CVMix/bin/cvmix < input.nl
if [ $? != 0 ]; then
  echo "Error in execution!"
  exit 1
fi

if [ "$OUTPUT" == "TRUE" ]; then
  if [ "`grep ltest1 input.nl | cut -d. -f 2`" != "false" ]; then
    if [ "${USE_NETCDF}" == "netcdf" ]; then
      ncdump -v zt,Ri_bulk test1.nc
    else
      cat test1.out
    fi
  fi

  if [ "`grep ltest4 input.nl | cut -d. -f 2`" != "false" ]; then
    if [ "${USE_NETCDF}" == "netcdf" ]; then
      ncdump -v zw,diff test4.nc
    else
      cat test4.out
    fi
  fi
fi
