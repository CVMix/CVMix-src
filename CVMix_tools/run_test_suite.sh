#!/bin/bash

# Convert shell return code to "PASS" or "FAIL"
# (0 = PASS, all other return codes = FAIL)
function check_return() {
  if [ $1 -eq 0 ]; then
    echo "PASS"
  else
    echo "FAIL"
  fi

}

#################################################

# Output test results
function print_status() {
  TEST_CNT=$((TEST_CNT+1))
  if [ "${STATUS}" == "FAIL" ]; then
    FAIL_CNT=$((FAIL_CNT+1))
  fi
  echo "${TEST_CNT}. $1: ${STATUS}"
}

#################################################

###############
# Global Vars #
###############

CVMIX_ROOT=`(cd ..; pwd -P)`
OUTFILE=${CVMIX_ROOT}/.test.out
TEST_CNT=0
FAIL_CNT=0
echo "Test Results:" > $OUTFILE

#########
# TESTS #
#########

# Check for --already-ran-setup
if [ "$1" == "--already-ran-setup" ]; then
  echo "Checking for existence of bld/.CVMix_env"
  if [ ! -f ${CVMIX_ROOT}/bld/.CVMix_env ]; then
    echo "ERROR: must run cvmix_setup prior to this script"
    exit 1
  fi
fi

# Code consistency check
cd ${CVMIX_ROOT}/CVMix_tools
echo "$ ./code_consistency.py"
./code_consistency.py
STATUS=$(check_return $?)
print_status "CodeConsistency.py" >> $OUTFILE

# Run pylint (if installed)
command -v pylint &> /dev/null
if [ $? -eq 0 ]; then
  cd ${CVMIX_ROOT}/CVMix_tools
  echo "$ pylint --rcfile=pylintrc code_consistency.py"
  pylint --rcfile=pylintrc code_consistency.py
  STATUS=$(check_return $?)
  print_status "pylint" >> $OUTFILE
fi

# Clean Fortran Code
cd ${CVMIX_ROOT}/src
echo "$ make allclean"
make allclean
STATUS=$(check_return $?)
print_status "make allclean" >> $OUTFILE

# Build libcvmix.a
cd ${CVMIX_ROOT}/src
echo "$ make lib"
make lib
STATUS=$(check_return $?)
print_status "make libcvmix.a" >> $OUTFILE

# Build stand-alone executable (only if library built successfully)
if [ "${STATUS}" == "PASS" ]; then
  echo "$ make"
  make
  STATUS=$(check_return $?)
  print_status "make cvmix.exe without netcdf" >> $OUTFILE
fi

# Only test Fortran executable if build was successful
if [ "${STATUS}" == "PASS" ]; then
  # Bryan-Lewis test
  cd ${CVMIX_ROOT}/reg_tests/Bryan-Lewis
  echo "$ ./Bryan-Lewis-test.sh"
  ./Bryan-Lewis-test.sh
  STATUS=$(check_return $?)
  print_status "Bryan Lewis" >> $OUTFILE

  # Double Diffusion test
  cd ${CVMIX_ROOT}/reg_tests/double_diff
  echo "$ ./double_diff-test.sh"
  ./double_diff-test.sh
  STATUS=$(check_return $?)
  print_status "Double Diffusion" >> $OUTFILE

  # Shear test
  cd ${CVMIX_ROOT}/reg_tests/shear
  echo "$ ./shear-test.sh"
  ./shear-test.sh
  STATUS=$(check_return $?)
  print_status "Shear" >> $OUTFILE

  # KPP test
  cd ${CVMIX_ROOT}/reg_tests/kpp
  echo "$ ./kpp-test.sh"
  ./kpp-test.sh
  STATUS=$(check_return $?)
  print_status "KPP.sh" >> $OUTFILE

  # Build stand-alone executable with netCDF
  cd ${CVMIX_ROOT}/src
  echo "$ make USE_NETCDF=TRUE"
  make USE_NETCDF=TRUE
  STATUS=$(check_return $?)
  print_status "make cvmix.exe with netcdf" >> $OUTFILE

  # Only test Fortran executable with netcdf if build was successful
  if [ "${STATUS}" == "PASS" ]; then
    # Tidal (Simmons) test
    cd ${CVMIX_ROOT}/reg_tests/tidal-Simmons
    echo "$ ./Simmons-test.sh -nc"
    ./Simmons-test.sh -nc
    STATUS=$(check_return $?)
    print_status "Tidal (requires netCDF)" >> $OUTFILE
  fi

fi



echo "----"
cat $OUTFILE
rm -f $OUTFILE
echo ""
echo "${TEST_CNT} tests were run, and $FAIL_CNT failed."
exit ${FAIL_CNT}

