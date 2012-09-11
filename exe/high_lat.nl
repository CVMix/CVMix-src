&cvmix_nml
mixing    ='BryanLewis'
nlev      = 30
ocn_depth = 5250
/
&BryanLewis_nml
! Low  latitude: 6.50e-5, 1.15e-4, 4.5e-3, 2500
! High latitude: 7.50e-5, 0.95e-4, 4.5e-3, 2500
! Orig BL paper: 8.00e-5, 1.05e-4, 4.5e-3, 2500
!
vdc1 = 7.50d-5
vdc2 = 0.95d-4
linv = 4.5d-3
dpth = 2500.0d0
/
