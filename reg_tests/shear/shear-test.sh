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
  ncdump -v Tdiff data_LMD.nc
  ncdump -v Mdiff data_PP1d.nc
else
  cat data_LMD.out
  cat data_PP1d.out
fi
