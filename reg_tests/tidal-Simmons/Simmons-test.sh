#!/bin/bash

##################
# Usage Subroutine
##################

usage() {
  echo './Simmons-test.sh [-nc|--netcdf]'
  echo ''
  echo 'By default, output is ascii'
  echo '  -h     | --help     Display this message'
  echo '  -nc    | --netcdf   Output to a netcdf file (data.nc) rather than ascii file'
  echo '  -clean | --clean    Delete all output data (data.nc and data.out files)'
}

#############
# Main Script
#############

while [ $# -gt 0 ]; do
  case $1 in
    -h|--help)
      usage
      exit 0
      ;;
    -nc|--netcdf)
      USE_NETCDF=netcdf
      ;;
    -clean|--clean)
      echo "rm -f single_col.nc diff.nc"
      rm -f single_col.nc diff.nc
      exit
      ;;
    * )
      echo "ERROR: improper use of script. Run ./Simmons-test.sh -h for usage."
      exit 1
      ;;
  esac
  shift
done

if [ "${USE_NETCDF}" != "netcdf" ]; then
  echo "Note: this test requires netCDF... use -nc flag"
  exit 1
fi

CVMix=$( cd ../.. ; pwd )
DATA_DIR=$CVMix/inputdata

# Files needed to run tidal example
GRID_FILE=gx1v6_130522.nc
PHYS_FILE=gx1v6_physics_130523.nc
ENERGY_FILE=tidal_energy_gx1v6_20130512.nc

# Input Data repository (currently using "svn export")
DATA_REPO=https://github.com/CVMix/CVMix-data/trunk/

echo "Checking for necessary input data files in ${DATA_DIR}..."
echo ""

echo "Looking for grid file ${GRID_FILE}..."
if [ -e ${DATA_DIR}/${GRID_FILE} ]; then
  echo "Found!"
else
  svn export ${DATA_REPO}/${GRID_FILE} ${DATA_DIR}
  echo "... Downloaded!"
fi

echo "Looking for physics file ${PHYS_FILE}..."
if [ -e ${DATA_DIR}/${PHYS_FILE} ]; then
  echo "Found!"
else
  svn export ${DATA_REPO}/${PHYS_FILE} ${DATA_DIR}
  echo "... Downloaded!"
fi

echo "Looking for energy map file ${ENERGY_FILE}..."
if [ -e ${DATA_DIR}/${ENERGY_FILE} ]; then
  echo "Found!"
else
  svn export ${DATA_REPO}/${ENERGY_FILE} ${DATA_DIR}
  echo "... Downloaded!"
fi

make -f $CVMix/src/Makefile CVMIX_ROOT=$CVMix ${USE_NETCDF}
if [ $? != 0 ]; then
  echo "Build error!"
  exit 1
fi

$CVMix/bin/cvmix < input.nl
if [ $? != 0 ]; then
  echo "Error in execution!"
  exit 1
fi

ncdump -v Tdiff single_col.nc
