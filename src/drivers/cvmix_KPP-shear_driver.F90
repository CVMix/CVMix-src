!BOP

! !ROUTINE: cvmix_KPPshear_driver

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

Program cvmix_KPPshear_driver

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                 &
                                    cvmix_data_type,          &
                                    cvmix_global_params_type, &
                                    cvmix_shear_params_type
  use cvmix_shear,           only : cvmix_init_shear,         &
                                    cvmix_coeffs_shear
  use cvmix_put_get,         only : cvmix_put
  use cvmix_output,          only : cvmix_output_open,        &
                                    cvmix_output_write,       &
                                    cvmix_output_close
!EOP
!BOC

  Implicit None

  type(cvmix_data_type)          :: CVmix_vars
  type(cvmix_global_params_type) :: CVmix_params
  type(cvmix_shear_params_type)  :: CVmix_KPP_params

  real(cvmix_r8), dimension(:),   allocatable, target :: Ri_g
  real(cvmix_r8), dimension(:),   allocatable, target :: viscosity
  real(cvmix_r8), dimension(:,:), allocatable, target :: diffusivity

  ! array index
  integer :: kw
  ! file index
  integer :: fid

  ! Namelist variables
  ! 1) General mixing parameters
  integer        :: nlev      ! number of Ri points to sample
  ! 2) KPP mixing parameters for column
  real(cvmix_r8) :: KPP_nu_zero, KPP_Ri_zero, KPP_exp

  ! Namelists that may be read in, depending on desired mixing scheme
  namelist/cvmix_nml/nlev
  namelist/KPP_nml/KPP_nu_zero, KPP_Ri_zero, KPP_exp

  ! Read general mixing parameters
  read(*, nml=cvmix_nml)

  ! Ri_g should increase from 0 to 1
  allocate(Ri_g(nlev+1))
  Ri_g(1) = 0.0_cvmix_r8
  do kw = 2,nlev+1
    Ri_g(kw) = Ri_g(kw-1) + 1.0_cvmix_r8/dble(nlev)
  end do

  ! Allocate memory to store viscosity and diffusivity values
  allocate(diffusivity(nlev+1,1), viscosity(nlev+1))

  ! Initialization for CVMix data types
  call cvmix_put(CVmix_params,  'max_nlev', nlev)
  call cvmix_put(CVmix_params,  'prandtl',  0.0_cvmix_r8)
  call cvmix_put(CVmix_vars,    'nlev',     nlev)
  ! Point CVmix_vars values to memory allocated above
  CVmix_vars%visc_iface => viscosity
  CVmix_vars%diff_iface => diffusivity
  CVmix_vars%Ri_iface => Ri_g

  ! Read / set KPP parameters
  read(*, nml=KPP_nml)
  call cvmix_init_shear(CVmix_KPP_params, 'KPP',                           &
                        KPP_nu_zero=KPP_nu_zero, KPP_Ri_zero=KPP_Ri_zero,  &
                        KPP_exp=KPP_exp)
  call cvmix_coeffs_shear(CVmix_vars, CVmix_KPP_params)

  ! Output
  ! data will have diffusivity from both columns (needed for NCL script)
#ifdef _NETCDF
  call cvmix_output_open(fid, "data.nc", "nc")
#else
  call cvmix_output_open(fid, "data.out", "ascii")
#endif

  call cvmix_output_write(fid, CVmix_vars, (/"Ri  ", "visc"/))

  call cvmix_output_close(fid)

!EOC

End program cvmix_KPPshear_driver
