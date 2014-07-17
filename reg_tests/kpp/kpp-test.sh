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
      ncdump -v zt,BulkRichardson test1.nc
    else
      echo "Test 1"
      echo "------"
      cat test1.out
    fi
  fi

  if [ "`grep ltest4 input.nl | cut -d. -f 2`" != "false" ]; then
    if [ "${USE_NETCDF}" == "netcdf" ]; then
      ncdump -v zw,Mdiff,Sdiff,Tdiff test4a.nc
      ncdump -v zw,Mdiff,Sdiff,Tdiff test4b.nc
    else
      echo "Test 4a"
      echo "------"
      cat test4a.out
      echo "Test 4b"
      echo "------"
      cat test4b.out
    fi
  fi

  if [ "`grep ltest5 input.nl | cut -d. -f 2`" != "false" ]; then
    if [ "${USE_NETCDF}" == "netcdf" ]; then
      ncdump -v zt,BulkRichardson test5.nc
    else
      echo "Test 5"
      echo "------"
      cat test5.out
    fi
  fi
fi
