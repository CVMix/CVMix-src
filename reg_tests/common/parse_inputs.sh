#/bin/bash

parse_inputs () {

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
      files=`ls *.out *.nc 2> /dev/null`
      if [ "$files" != "" ]; then
        clean_cmd="rm -f $files"
        echo $clean_cmd
      else
        clean_cmd='echo No output files to remove!'
      fi
      $clean_cmd
      exit
    ;;
    -nb|--nobuild)
      if [ ! -e $CVMix/bin/cvmix ]; then
        echo "ERROR: executable not found!"
        exit 1
      fi
      NO_BUILD=1
    ;;
    * )
      if [ "$ThisScript" == "BL_test.sh" ]; then
        if [ "$1" == "-mc" ] || [ "$1" == "-memcopy" ]; then
          DRIVER=mem_copy
        else
          echo "ERROR: $1 is not a valid argument for $ThisScript." \
               "Run '$ThisScript -h' for proper usage."
          exit 1
        fi
      else
        echo "ERROR: $1 is not a valid argument for $ThisScript." \
             "Run '$ThisScript -h' for proper usage."
        exit 1
      fi
    ;;
  esac
  shift
done

}
