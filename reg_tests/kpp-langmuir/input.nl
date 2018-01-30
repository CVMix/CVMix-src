&cvmix_nml
mix_type = 'kpp-langmuir'
/
 
! Langmuir Test params
&langmuir_col_nml
  nlev           = 256
  max_nlev       = 256
  layer_thick    = 0.64d0
  interp_type    = 'linear'
  ri_crit        = 0.30d0
  llangmuir_efactor = .false.
  lamult         = 1.00d0
  b0             = -2.3357e-09
  b0sol          = 0.0000e+00
  jerlov_water_type = 3
  ustar          = 6.0000e-03
  Coriolis       = 1.0284e-04
  infile         = '/Users/qingli/GoogleDrive/BrownWork/matlab/ncarles/testCVMix/uvbPfl.txt'
  outfile        = 'langmuir_off.out'
  lnoDGat1       = .false.
/
