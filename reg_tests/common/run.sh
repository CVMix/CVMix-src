#!/bin/bash

run () {

  if [ "${CMAKE_BUILD}" == "TRUE" ]; then
    CVMIX_EXE=$CVMix/bld/cmake_bld/cvmix_driver
  else
    CVMIX_EXE=$CVMix/bin/cvmix
  fi
  ${CVMIX_EXE} < $NAMELIST

  if [ $? != 0 ]; then
    echo "Error in execution!"
    exit 3
  fi

}
