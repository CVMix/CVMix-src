#!/bin/bash

check_inputdata () {

  DATA_DIR=$CVMix/inputdata
  DATA_REPO=https://github.com/CVMix/CVMix-data/trunk/

  for file in "$@"; do
    echo "Looking for $file..."
    if [ -e ${DATA_DIR}/$file ]; then
      echo "Found!"
    else
      svn export ${DATA_REPO}/$file ${DATA_DIR}/$file
      echo "... Downloaded!"
    fi
  done

}
