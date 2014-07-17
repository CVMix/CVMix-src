#!/bin/bash

check_inputdata () {

  DATA_REPO=https://github.com/CVMix/CVMix-data/trunk/

  for file in "$@"; do
    echo "Looking for $file..."
    if [ -e ${INPUTDATA_DIR}/$file ]; then
      echo "Found!"
    else
      svn export ${DATA_REPO}/$file ${INPUTDATA_DIR}/$file
      echo "... Downloaded!"
    fi
  done

}
