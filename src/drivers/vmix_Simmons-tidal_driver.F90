!BOP

! !ROUTINE: vmix_Simmonstidal_driver

! !DESCRIPTION: A stand-alone driver for the CVMix package. This particular
!  driver generates the tidal-mixing coefficient defined in Equation (??) of
!  Simmons, et al., in a single column and then outputs the column to allow
!  recreation of Figure ? from the same paper.
!\\
!\\

! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !INTERFACE:

Program vmix_Simmonstidal_driver

! !USES:

  use vmix_kinds_and_types, only : vmix_r8,                 &
                                   vmix_data_type,          &
                                   vmix_global_params_type, &
                                   vmix_tidal_params_type
  use vmix_tidal,           only : vmix_init_tidal,         &
                                   vmix_coeffs_tidal
  use vmix_put_get,         only : vmix_put
  use vmix_output,          only : vmix_output_open,        &
                                   vmix_output_write,       &
                                   vmix_output_close
!EOP
!BOC

  Implicit None

  type (vmix_data_type)          :: Vmix_vars
  type (vmix_global_params_type) :: Vmix_params
  type (vmix_tidal_params_type)  :: Vmix_Simmons_params

  real(kind=vmix_r8), dimension(:),   allocatable, target :: viscosity
  real(kind=vmix_r8), dimension(:,:), allocatable, target :: diffusivity

  ! file index
  integer :: fid

  ! Namelist variables
  ! 1) General mixing parameters
  integer                    :: nlev      ! number of Ri points to sample
  ! 2) Other variables for tidal mixing go here

  ! Namelists that may be read in, depending on desired mixing scheme
  namelist/cvmix_nml/nlev
  ! namelist/Simmons_nml/

  ! Read general mixing parameters
  read(*, nml=cvmix_nml)

  ! Allocate memory to store viscosity and diffusivity values
  allocate(diffusivity(nlev+1,1), viscosity(nlev+1))

  ! Initialization for CVMix data types
  call vmix_put(Vmix_params,  'max_nlev', nlev)
  call vmix_put(Vmix_params,  'prandtl',  0.0d0)
  call vmix_put(Vmix_vars, 'nlev', nlev)
  ! Point Vmix_vars values to memory allocated above
  Vmix_vars%visc_iface => viscosity
  Vmix_vars%diff_iface => diffusivity

  ! Read / set Simmons parameters
!  read(*, nml=KPP_nml)
  call vmix_init_tidal(Vmix_Simmons_params, 'simmons')
  call vmix_coeffs_tidal(Vmix_vars, Vmix_Simmons_params)

  ! Output
  ! data will have diffusivity from both columns (needed for NCL script)
#ifdef _NETCDF
  call vmix_output_open(fid, "data.nc", "nc")
#else
  call vmix_output_open(fid, "data.out", "ascii")
#endif

  call vmix_output_write(fid, Vmix_vars, (/"diff", "visc"/))

  call vmix_output_close(fid)

!EOC

End program vmix_Simmonstidal_driver
