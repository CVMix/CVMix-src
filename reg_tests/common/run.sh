#!/bin/bash

run () {

  if [ $ThisScript == "BL_test.sh" ]; then
    if [ "$DRIVER" == "mem_copy" ]; then
      $CVMix/bin/cvmix < input_memcopy.nl
    else
      $CVMix/bin/cvmix < input_pointer.nl
    fi
  else
    $CVMix/bin/cvmix < input.nl
  fi

  if [ $? != 0 ]; then
    echo "Error in execution!"
    exit 3
  fi

}
