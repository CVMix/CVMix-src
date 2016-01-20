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

# If different inputdata directory specified, update namelist
if [ "$INPUTDATA_DIR" != "$CVMix/inputdata" ]; then
  rm -f input2.nl
  NAMELIST=input2.nl
  sed s+grid_file.\*+"grid_file = \'$INPUTDATA_DIR/$GRID_FILE\'"+ ./input.nl      |\
  sed s+physics_file.\*+"physics_file = \'$INPUTDATA_DIR/$PHYS_FILE\'"+           |\
  sed s+energy_flux_file.\*+"energy_flux_file = \'$INPUTDATA_DIR/$ENERGY_FILE\'"+  \
      > input2.nl
fi

build
run

# (2) Look at output
ncdump -v Mdiff,Tdiff single_col.nc
