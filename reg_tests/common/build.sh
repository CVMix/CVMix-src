#!/bin/bash

build () {
  if [ -z ${NO_BUILD} ]; then
    make -f $CVMix/src/Makefile CVMIX_ROOT=$CVMix ${USE_NETCDF}
    if [ $? != 0 ]; then
      echo "Build error!"
      exit 2
    fi
  else
    if [ -e ../../bld/.netcdf_info ]; then
      if [ "`cat ../../bld/.netcdf_info`" == "YES" ]; then
        USE_NETCDF=netcdf
      fi
    fi  
  fi

}
