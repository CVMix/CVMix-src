!BOP
!\newpage
! !ROUTINE: cvmix_tidal_driver

! !DESCRIPTION: A routine to test the Simmons implementation of tidal mixing.
!\\
!\\

! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !INTERFACE:

Subroutine cvmix_tidal_driver(nlev)

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                 &
                                    cvmix_strlen,             &
                                    cvmix_data_type,          &
                                    cvmix_global_params_type, &
                                    cvmix_tidal_params_type
  use cvmix_tidal,           only : cvmix_init_tidal,         &
                                    cvmix_coeffs_tidal
  use cvmix_put_get,         only : cvmix_put
  use cvmix_io,              only : cvmix_io_open,            &
                                    cvmix_input_read,         &
                                    cvmix_output_write,       &
                                    cvmix_io_close

  Implicit None

! !INPUT PARAMETERS
  integer, intent(in) :: nlev

!EOP
!BOC

  ! CVMix datatypes
  type(cvmix_data_type)          :: CVmix_vars
  type(cvmix_global_params_type) :: CVmix_params
  type(cvmix_tidal_params_type)  :: CVmix_Simmons_params

  real(cvmix_r8), dimension(:),   allocatable, target :: viscosity
  real(cvmix_r8), dimension(:,:), allocatable, target :: diffusivity

  ! file index
  integer :: fid

  ! Namelist variables
  character(len=cvmix_strlen) :: grid_file, energy_flux_file, energy_flux_var
  integer :: nlon, nlat

  ! Local variables
  real(cvmix_r8), dimension(:,:), allocatable :: ocn_depth, energy_flux
  real(cvmix_r8) :: depth_fill, flux_fill, my_min, my_max
  integer :: i,j

  ! Namelists that may be read in, depending on desired mixing scheme
  namelist/Simmons_nml/grid_file, energy_flux_file, energy_flux_var, nlon, nlat

  ! Allocate memory to store viscosity and diffusivity values
  allocate(diffusivity(nlev+1,1), viscosity(nlev+1))

  ! Initialization for CVMix data types
  call cvmix_put(CVmix_params, 'max_nlev',          nlev)
  call cvmix_put(CVmix_params,  'prandtl',  0.0_cvmix_r8)
  call cvmix_put(CVmix_vars,       'nlev',          nlev)
  call cvmix_put(CVmix_vars,   'surf_hgt',  0.0_cvmix_r8)
  ! Point CVmix_vars values to memory allocated above
  CVmix_vars%visc_iface => viscosity
  CVmix_vars%diff_iface => diffusivity

  ! Read / set Simmons parameters
  ! Default values
  grid_file = "none"
  energy_flux_file = "none"
  energy_flux_var = "none"
  nlon = 320
  nlat = 184
  read(*, nml=Simmons_nml)
  allocate(energy_flux(nlon, nlat), ocn_depth(nlon, nlat))
  call cvmix_io_open(fid, trim(grid_file), 'nc', read_only=.true.)
  call cvmix_input_read(fid, 'H', ocn_depth)
  call cvmix_io_close(fid)
  call cvmix_io_open(fid, trim(energy_flux_file), 'nc', read_only=.true.)
  call cvmix_input_read(fid, trim(energy_flux_var), energy_flux)
  call cvmix_io_close(fid)
  ! Note: at this time, not ignoring missing value
  print*, "Min and Max of ocean depth:"
  depth_fill = maxval(ocn_depth)
  my_min = 99999.0_cvmix_r8
  my_max = 0.0_cvmix_r8
  do i=1,nlon
    do j=1,nlat
      if (ocn_depth(i,j).ne.depth_fill) then
        my_min = min(ocn_depth(i,j),my_min)
        my_max = max(ocn_depth(i,j),my_max)
      end if
    end do
  end do
  print*, my_min, my_max

  print*, "Min and Max of energy flux:"
  flux_fill = minval(energy_flux)
  my_min = 99999.0_cvmix_r8
  my_max = 0.0_cvmix_r8
  do i=1,nlon
    do j=1,nlat
      if (energy_flux(i,j).ne.flux_fill) then
        my_min = min(energy_flux(i,j),my_min)
        my_max = max(energy_flux(i,j),my_max)
      end if
    end do
  end do
  print*, my_min, my_max

  call cvmix_init_tidal(CVmix_Simmons_params, 'Simmons', 'mks')
  print*, "Namelist variables:"
  print*, "mix_scheme = ", trim(CVmix_Simmons_params%mix_scheme)
  print*, "efficiency = ", CVmix_Simmons_params%efficiency
  print*, "vertical_decay_scale = ", CVmix_Simmons_params%vertical_decay_scale
  print*, "max_coefficient = ", CVmix_Simmons_params%max_coefficient
  print*, "local_mixing_frac = ", CVmix_Simmons_params%local_mixing_frac
  print*, "depth_cutoff = ", CVmix_Simmons_params%depth_cutoff
  ! Picking a random column to test (for now)
  print*, "Depth: ", ocn_depth(228, 125)
  print*, "Flux: ", energy_flux(228, 125)
  CVmix_vars%ocn_depth = ocn_depth(228, 125)
  call cvmix_coeffs_tidal(CVmix_vars, CVmix_Simmons_params)

  ! Output
  ! data will have diffusivity from both columns (needed for NCL script)
#ifdef _NETCDF
  call cvmix_io_open(fid, "data.nc", "nc")
#else
  call cvmix_io_open(fid, "data.out", "ascii")
#endif

  call cvmix_output_write(fid, CVmix_vars, (/"diff", "visc"/))

  call cvmix_io_close(fid)

!EOC

End Subroutine cvmix_tidal_driver
