#!/bin/bash

# (1) Load required routines
NAMELIST=input_pointer.nl
. ../common/environ.sh
. ../common/usage.sh
. ../common/parse_inputs.sh
. ../common/build.sh
. ../common/run.sh

DRIVER=pointers

parse_inputs $@
build
run

# (2) Look at output
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
