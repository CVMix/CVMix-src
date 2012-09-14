#!/bin/bash

CVMix=$PWD/../..
make -f $CVMix/bld/Makefile VMIX_DRIVER=vmix_BL_driver-mem_copy.F90 VMIX_ROOT=$CVMix

$CVMix/bld/exe/cvmix < input.nl
