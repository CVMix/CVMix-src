#!/bin/bash

# (1) Load required routines
. ../common/environ.sh
. ../common/usage.sh
. ../common/parse_inputs.sh
. ../common/build.sh
. ../common/run.sh
. ../common/check_inputdata.sh

parse_inputs $@
if [ "${USE_NETCDF}" != "netcdf" ]; then
  echo "Note: this test requires netCDF... use -nc flag"
  exit 1
fi

# Files needed to run tidal example
GRID_FILE=gx1v6_130522.nc
PHYS_FILE=gx1v6_physics_130523.nc
ENERGY_FILE=tidal_energy_gx1v6_20130512.nc

# Input Data repository (currently using "svn export")
check_inputdata $GRID_FILE $PHYS_FILE $ENERGY_FILE

build
run

# (2) Look at output
ncdump -v Tdiff single_col.nc
