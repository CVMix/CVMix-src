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
                                    cvmix_data_type,          &
                                    cvmix_global_params_type, &
                                    cvmix_tidal_params_type
  use cvmix_tidal,           only : cvmix_init_tidal,         &
                                    cvmix_coeffs_tidal
  use cvmix_put_get,         only : cvmix_put
  use cvmix_io,              only : cvmix_output_open,        &
                                    cvmix_output_write,       &
                                    cvmix_output_close

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
  ! Variables for tidal mixing go here

  ! Namelists that may be read in, depending on desired mixing scheme
  ! namelist/Simmons_nml/

  ! Allocate memory to store viscosity and diffusivity values
  allocate(diffusivity(nlev+1,1), viscosity(nlev+1))

  ! Initialization for CVMix data types
  call cvmix_put(CVmix_params,  'max_nlev', nlev)
  call cvmix_put(CVmix_params,  'prandtl',  0.0_cvmix_r8)
  call cvmix_put(CVmix_vars,    'nlev',     nlev)
  ! Point CVmix_vars values to memory allocated above
  CVmix_vars%visc_iface => viscosity
  CVmix_vars%diff_iface => diffusivity

  ! Read / set Simmons parameters
!  read(*, nml=Simmons_nml)
  call cvmix_init_tidal(CVmix_Simmons_params, 'simmons')
  call cvmix_coeffs_tidal(CVmix_vars, CVmix_Simmons_params)

  ! Output
  ! data will have diffusivity from both columns (needed for NCL script)
#ifdef _NETCDF
  call cvmix_output_open(fid, "data.nc", "nc")
#else
  call cvmix_output_open(fid, "data.out", "ascii")
#endif

  call cvmix_output_write(fid, CVmix_vars, (/"diff", "visc"/))

  call cvmix_output_close(fid)

!EOC

End Subroutine cvmix_tidal_driver
