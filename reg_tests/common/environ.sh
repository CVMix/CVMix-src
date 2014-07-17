#!/bin/bash

CVMix=$( cd ../.. ; pwd )
ThisScript=$(basename $0)
if [ -z $NAMELIST ]; then
  NAMELIST=input.nl
fi
INPUTDATA_DIR=$CVMix/inputdata
