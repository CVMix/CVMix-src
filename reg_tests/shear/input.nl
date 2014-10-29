&cvmix_nml
mix_type = 'shear'
nlev     = 30
max_nlev = 60
/
! Shear mixing parameters from LMD94 paper
&LMD_nml
LMD_nu_zero = 5d-3
LMD_Ri_zero = 0.7d0
LMD_exp     = 3.0d0
/
! Shear mixing parameters from PP81 paper
&PP_nml
PP_nu_zero = 5d-3
PP_alpha   = 5.0d0
PP_exp     = 2.0d0
/
