#!/bin/bash

build () {
  if [ -z ${NO_BUILD} ]; then
    if [ "${CMAKE_BUILD}" == "TRUE" ]; then
      CWD=$PWD
      cd $CVMix/bld/cmake_bld
      cmake $CVMix -DCVMIX_BUILD_DRIVER=on \
                   -DCMAKE_Fortran_COMPILER=gfortran \
                   -DCMAKE_INSTALL_PREFIX=$CVMix/cmake_bin \
      && make && make install
      STATUS=$?
      cd $CWD
    else
      make -f $CVMix/src/Makefile CVMIX_ROOT=$CVMix ${USE_NETCDF}
      STATUS=$?
    fi
    if [ $STATUS != 0 ]; then
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
