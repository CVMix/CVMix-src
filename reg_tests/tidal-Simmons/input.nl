&cvmix_nml
mix_type = 'tidal'
nlev    = 30
/
! Tidal mixing parameters from Simmons paper
&Simmons_nml
grid_file = '../../inputdata/gx1v6_130522.nc'
physics_file = '../../inputdata/gx1v6_physics_130523.nc'
energy_flux_file = '../../inputdata/tidal_energy_gx1v6_20130512.nc'
energy_flux_var = 'TIDAL_ENERGY_FLUX'
nlon = 320
nlat = 384
/
