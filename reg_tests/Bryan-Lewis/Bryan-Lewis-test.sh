#!/bin/bash

# (1) Load required routines
. ../common/environ.sh
. ../common/usage.sh
. ../common/parse_inputs.sh
. ../common/build.sh
. ../common/run.sh

parse_inputs $@
build
run

# (2) Look at output
if [ "${USE_NETCDF}" == "netcdf" ]; then
  if [ -f data_pointer.nc ] && [ -f data_memcopy.nc ]; then
    diff data_pointer.nc data_memcopy.nc
    if [ $? -ne 0 ]; then
      echo "WARNING: output files are different, they should be identical!"
      ncdump -v Tdiff data_pointer.nc
      ncdump -v Tdiff data_memcopy.nc
      echo "WARNING: output files are different, they should be identical!"
    else
      ncdump -v Tdiff data_pointer.nc
    fi
  else
    echo "ERROR: Output not found!"
    exit 4
  fi
else
  if [ -f data_pointer.out ] && [ -f data_memcopy.out ]; then
    diff data_pointer.out data_memcopy.out
    if [ $? -ne 0 ]; then
      echo "WARNING: output files are different, they should be identical!"
      cat data_pointer.out
      echo "---"
      cat data_memcopy.out
      echo "WARNING: output files are different, they should be identical!"
    else
      cat data_pointer.out
    fi
  else
    echo "ERROR: Output not found!"
    exit 4
  fi
fi
