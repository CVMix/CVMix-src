#!/bin/bash -e
# Since github dropped support for svn hooks, now need to set up
# a sparse checkout to get inputdata.
# NOTE: this script must be run from ${CVMIX_ROOT}!

cd inputdata
git init
git config core.sparseCheckout true
git remote add -f CVMix-data https://github.com/CVMix/CVMix-data.git