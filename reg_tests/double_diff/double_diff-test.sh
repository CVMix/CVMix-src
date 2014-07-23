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
  ncdump -v Tdiff,Sdiff data.nc
else
  cat data.out
fi
