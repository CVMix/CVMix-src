&cvmix_nml
mix_type = 'tidal'
/
! Tidal mixing parameters from Simmons paper
&Simmons_nml
grid_file = '../../inputdata/gx1v6_130522.nc'
physics_file = '../../inputdata/gx1v6_physics_130523.nc'
energy_flux_file = '../../inputdata/tidal_energy_gx1v6_20130512.nc'
energy_flux_var = 'TIDAL_ENERGY_FLUX'
! Aleutians (170.24 W, 51.04 N)
! lon_out = 203
! lat_out = 309
! Scotland (6.04 W, 58.85 N)
! lon_out = 35
! lat_out = 345
! Madagascar (49.44 E, 11.62 S)
lon_out = 80 
lat_out = 144
/
