&cvmix_nml
mix_type = 'kpp'
/

! Test 1 params
&kpp_col1_nml
!  ltest1         = .true.
  nlev1          = 4
  layer_thick    = 10.0d0
  interp_type_t1 = 'quadratic'
  hmix           = -15.0d0
  ri_crit        = 0.3d0
/

! Test 2 params
&kpp_col2_nml
!  ltest2 = .true.
/

! Test 3 params
&kpp_col3_nml
!  ltest3 = .true.
/

! Test 4 params
&kpp_col4_nml
  ltest4         = .true.
  interp_type_t4 = "quadratic"
  OBL_depth      = 14.0d0
  lnoDGat1       = .true.
/
