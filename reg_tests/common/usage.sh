#/bin/bash

usage() {
  if [ $ThisScript == "Simmons-test.sh" ]; then
    USE="$ThisScript -nc|--netcdf"
  else
    USE="$ThisScript [-nc|--netcdf]"
  fi
  USE="$USE [-nb|--no_build]"
  echo $USE
  echo ''
  if [ $ThisScript == "Simmons-test.sh" ]; then
    echo 'This test requires netcdf input, so you must run with the -nc option.'
  else
    echo 'By default, output is ascii.'
  fi
  echo '  -h     | --help      Display this message'
  echo '  -nc    | --netcdf    Output in netcdf format rather than ascii'
  echo '  -nb    | --no_build  Do not rebuild the executable'
  echo '  -clean | --clean     Delete all output data (*.nc and *.out files)'
}

