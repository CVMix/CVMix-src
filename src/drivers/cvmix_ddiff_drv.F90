!BOP
!\newpage
! !ROUTINE: cvmix_ddiff_driver

! !DESCRIPTION: A routine to test the double diffusion mixing module.
!\\
!\\

! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !INTERFACE:

Subroutine cvmix_ddiff_driver(nlev)

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                 &
                                    cvmix_data_type,          &
                                    cvmix_global_params_type, &
                                    cvmix_ddiff_params_type
  use cvmix_ddiff,           only : cvmix_init_ddiff,         &
                                    cvmix_coeffs_ddiff
  use cvmix_put_get,         only : cvmix_put
  use cvmix_output,          only : cvmix_output_open,        &
                                    cvmix_output_write,       &
                                    cvmix_output_close

  Implicit None

! !INPUT PARAMETERS:
  integer, intent(in) :: nlev

!EOP
!BOC

  ! CVMix datatypes
  type(cvmix_data_type)          :: CVmix_vars
  type(cvmix_global_params_type) :: CVmix_params
  type(cvmix_ddiff_params_type)  :: CVmix_ddiff_params

  real(cvmix_r8), dimension(:),   allocatable, target :: viscosity
  real(cvmix_r8), dimension(:,:), allocatable, target :: diffusivity

  ! file index
  integer :: fid

  ! Namelist variables
  ! Variables for double diffusion mixing go here

  ! Namelists that may be read in, depending on desired mixing scheme
  ! namelist/ddiff_nml/

  ! Allocate memory to store viscosity and diffusivity values
  allocate(diffusivity(nlev+1,1), viscosity(nlev+1))

  ! Initialization for CVMix data types
  call cvmix_put(CVmix_params,  'max_nlev', nlev)
  call cvmix_put(CVmix_params,  'prandtl',  0.0_cvmix_r8)
  call cvmix_put(CVmix_vars,    'nlev',     nlev)
  ! Point CVmix_vars values to memory allocated above
  CVmix_vars%visc_iface => viscosity
  CVmix_vars%diff_iface => diffusivity

  ! Read / set double diffusion parameters
!  read(*, nml=ddiff_nml)
  call cvmix_init_ddiff(CVmix_ddiff_params)
  call cvmix_coeffs_ddiff(CVmix_vars, CVmix_ddiff_params)

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

End Subroutine cvmix_ddiff_driver
