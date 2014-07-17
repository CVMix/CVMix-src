#!/bin/bash

run () {

  $CVMix/bin/cvmix < $NAMELIST

  if [ $? != 0 ]; then
    echo "Error in execution!"
    exit 3
  fi

}
