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
  interp_type_t4 = "quadratic"
  OBL_depth      = 14.0d0
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
  lavg_N_or_Nsqr = .true.
/
