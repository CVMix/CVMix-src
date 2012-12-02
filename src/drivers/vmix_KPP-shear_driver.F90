!BOP

! !ROUTINE: vmix_KPPshear_driver

! !DESCRIPTION: A stand-alone driver for the CVMix package. This particular
!  driver generates the shear-mixing coefficient defined in Equation (28) of
!  Large, et al., in a single column and then outputs the column to allow
!  recreation of Figure 3 from the same paper. Note that here each "level" of
!  the column denotes a different local gradient Richardson number rather than
!  a physical ocean level. All memory is declared in the driver, and the CVMix
!  data type points to the local variables.
!\\
!\\

! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !INTERFACE:

Program vmix_KPPshear_driver

! !USES:

  use vmix_kinds_and_types
  use vmix_shear
  use vmix_put_get
  use vmix_output
!EOP
!BOC

  Implicit None

  type (vmix_data_type)          :: Vmix_vars
  type (vmix_global_params_type) :: Vmix_params
  type (vmix_shear_params_type)  :: Vmix_KPP_params

  real(kind=vmix_r8), dimension(:),   allocatable, target :: Ri_g
  real(kind=vmix_r8), dimension(:),   allocatable, target :: viscosity
  real(kind=vmix_r8), dimension(:,:), allocatable, target :: diffusivity

  ! array index
  integer :: kw
  ! file index
  integer :: fid

  ! Namelist variables
  ! 1) General mixing parameters
  integer                    :: nlev      ! number of Ri points to sample
  ! 2) KPP mixing parameters for column
  real(kind=vmix_r8)         :: nu_zero, Ri_zero, p_one

  ! Namelists that may be read in, depending on desired mixing scheme
  namelist/cvmix_nml/nlev
  namelist/KPP_nml/nu_zero, Ri_zero, p_one

  ! Read general mixing parameters
  read(*, nml=cvmix_nml)

  ! Ri_g should increase from 0 to 1
  allocate(Ri_g(nlev+1))
  Ri_g(1) = 0.0d0
  do kw = 2,nlev+1
    Ri_g(kw) = Ri_g(kw-1) + 1.0d0/dble(nlev)
  end do

  ! Allocate memory to store viscosity and diffusivity values
  allocate(diffusivity(nlev+1,1), viscosity(nlev+1))

  ! Initialization for CVMix data types
  call vmix_put(Vmix_params,  'max_nlev', nlev)
  call vmix_put(Vmix_params,  'prandtl',  0.0d0)
  call vmix_put(Vmix_vars, 'nlev', nlev)
  ! Point Vmix_vars values to memory allocated above
  Vmix_vars%visc_iface => viscosity
  Vmix_vars%diff_iface => diffusivity
  Vmix_vars%Ri_t_iface => Ri_g
  Vmix_vars%z_iface    => Ri_g ! For output purposes!

  ! Read / set KPP parameters
  read(*, nml=KPP_nml)
  call vmix_init_shear(Vmix_KPP_params, 'KPP',                              &
                       nu_zero = nu_zero, Ri_zero = Ri_zero, p_one = p_one)
  call vmix_coeffs_shear(Vmix_vars, Vmix_KPP_params)

  ! Output
  ! data will have diffusivity from both columns (needed for NCL script)
#ifdef _NETCDF
  call vmix_output_open(fid, "data.nc", "nc")
#else
  call vmix_output_open(fid, "data.out", "ascii")
#endif

  call vmix_output_write(fid, Vmix_vars, (/"Ri  ", "visc"/))

  call vmix_output_close(fid)

!EOC

End program vmix_KPPshear_driver
