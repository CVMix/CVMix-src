&cvmix_nml
mix_type  = 'BryanLewis'
nlev      = 30
max_nlev  = 60
ocn_depth = 5250.0d0
/
! BL Parameters --
! Low  latitude: 6.50e-5, 1.15e-4, 4.5e-3, 2500
! High latitude: 7.50e-5, 0.95e-4, 4.5e-3, 2500
! Orig BL paper: 8.00e-5, 1.05e-4, 4.5e-3, 2500
!
&BryanLewis1_nml
col1_vdc1 = 6.50d-5
col1_vdc2 = 1.15d-4
col1_linv = 4.5d-3
col1_dpth = 2500.0d0
/
&BryanLewis2_nml
col2_vdc1 = 7.50d-5
col2_vdc2 = 0.95d-4
col2_linv = 4.5d-3
col2_dpth = 2500.0d0
/
