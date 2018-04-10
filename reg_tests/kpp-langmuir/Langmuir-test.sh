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
if [ "$OUTPUT" == "TRUE" ]; then
  if [ "`grep ltest1 input.nl | cut -d. -f 2`" != "false" ]; then
    if [ "${USE_NETCDF}" == "netcdf" ]; then
      ncdump -v zt,Ri_bulk langmuir_test.nc
    else
      echo "Langmuir Test"
      echo "------"
      cat langmuir_test.out
    fi
  fi
fi
