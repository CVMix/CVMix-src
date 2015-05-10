&cvmix_nml
mix_type = 'kpp'
/

! Test 1 params
&kpp_col1_nml
  ltest1         = .true.
  nlev1          = 4
  layer_thick1   = 10.0d0
  interp_type_t1 = 'quadratic'
  hmix1          = -15.0d0
  ri_crit        = 0.3d0
/

! Test 2 params
&kpp_col2_nml
  ltest2 = .true.
/

! Test 3 params
&kpp_col3_nml
  ltest3 = .true.
/

! Test 4 params
&kpp_col4_nml
  ltest4         = .true.
  interp_type_t4 = 'quadratic'
  OBL_levid4     = 3
  lnoDGat1       = .true.
/

! Test 5 params
&kpp_col5_nml
  ltest5         = .true.
  nlev5          = 10
  layer_thick5   = 5.0d0
  hmix5          = 17.0d0
  ! Parameter settings to match LMD94 (linear interp, average Nsqr)
  interp_type_t5 = "linear"
/

! Test 6 params
&kpp_col6_nml
  ltest6     = .true.
  vonkarman6 = 0.4d0
  tao        = 0.2d0
  rho0       = 1035.0d0
  grav       = 9.8d0
  alpha      = 2.5d-4
  Qnonpen    = -100.0d0
  Cp0        = 3992.0d0
  OBL_depth6 = 6000.0d0
/
