#!/bin/bash -e

check_inputdata () {

  DATA_REPO=https://github.com/CVMix/CVMix-data/trunk/

  cd ${INPUTDATA_DIR}
  ALL_FOUND=TRUE
  for file in "$@"; do
    echo "Looking for $file..."
    if [ -e $file ]; then
      echo "Found!"
    else
      echo "${file}" >> .git/info/sparse-checkout
      ALL_FOUND=FALSE
      echo "... added to sparse-checkout!"
    fi
  done

  if [ "${ALL_FOUND}" == "FALSE" ]; then
    git checkout master
  fi

}