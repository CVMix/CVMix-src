#/bin/bash

usage() {
  if [ $ThisScript == "Simmons-test.sh" ]; then
    USE="$ThisScript -nc|--netcdf"
  else
    USE="$ThisScript [-nc|--netcdf]"
  fi
    USE="$USE [-nb|--no_build]"
  if [ $ThisScript == "BL_test.sh" ]; then
    USE="$USE [-mc|--memcopy]"
  fi
  echo $USE
  echo ''
  if [ $ThisScript == "BL_test.sh" ]; then
    echo 'By default, this script allocates memory in the driver and the CVMix datatypes point back to memory. Output is ascii.'
  else
  if [ $ThisScript == "Simmons-test.sh" ]; then
    echo 'This test requires netcdf input, so you must run with the -nc option.'
  else
    echo 'By default, output is ascii.'
  fi
  fi
  echo '  -h     | --help      Display this message'
  if [ $ThisScript == "BL_test.sh" ]; then
    echo '  -mc    | --memcopy   Allocate memory in the CVMix datatypes rather than point'
  fi
  echo '  -nc    | --netcdf    Output in netcdf format rather than ascii'
  echo '  -nb    | --no_build  Do not rebuild the executable'
  echo '  -clean | --clean     Delete all output data (*.nc and *.out files)'
}

