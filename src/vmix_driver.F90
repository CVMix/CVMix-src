!BOI

! !TITLE: In-code documentation for CVMix
! !AUTHORS: Many contributors from GFDL, LANL, and NCAR
! !AFFILIATION: GFDL, LANL, and NCAR
! !DATE: \today

!EOI
!BOP

! !ROUTINE: vmix_driver

! !DESCRIPTION: The stand-alone driver for the CVMix package. The type of
!  mixing comes from the cvmix\_nml namelist, and then other namelists
!  contain parameters for the specific mixing methods.
!\\
!\\

! !INTERFACE:

Program vmix_driver

! !USES:

  use vmix_kinds_and_types
  use vmix_background
  use vmix_convection
  use vmix_put_get
  use vmix_output

!EOP
!BOC

  Implicit None

  type (vmix_data_type)                         :: Vmix_vars
  type (vmix_global_params_type)                :: Vmix_params
  type (vmix_bkgnd_params_type)                 :: Vmix_BL_params
  real(kind=vmix_r8), dimension(:), allocatable, target :: iface_depth, &! nlev+1
                                                           viscosity     ! nlev+1
  real(kind=vmix_r8), dimension(:,:), allocatable, target :: diffusivity ! nlev+1 x 1

  ! array indices
  integer :: kw

  ! Namelist variables
  ! 1) General mixing parameters
  character(len=vmix_strlen) :: mixing    ! type of mixing
  integer                    :: nlev      ! number of levels for column
  real(kind=vmix_r8)         :: ocn_depth ! Depth of ocn
  ! 2) Bryan-Lewis mixing parameters
  ! k = vdc1 + (vdc2/PI)*atan((z-dpth)*linv)
  real(kind=vmix_r8)         :: vdc1, vdc2, linv, dpth

  ! Namelists that may be read in, depending on desired mixing scheme
  namelist/cvmix_nml/mixing, nlev, ocn_depth
  namelist/BryanLewis_nml/vdc1, vdc2, linv, dpth

  ! Read general mixing parameters
  read(*, nml=cvmix_nml)

  ! Calculate depth of cell interfaces based on number of levels and ocean
  ! depth (also allocate memory for diffusivity and viscosity)
  allocate(iface_depth(nlev+1), diffusivity(nlev+1,1), viscosity(nlev+1))
  iface_depth(1) = 0.0d0
  do kw = 2,nlev+1
    iface_depth(kw) = iface_depth(kw-1) + ocn_depth/dble(nlev)
  end do

  ! Initialization for CVMix data types
  ! Note that vmix_put will allocate memory for arrays as needed
  ! Alternatively, for visc, diff, and zw, you could point directly
  ! at the data
  call vmix_put(Vmix_params,  'max_nlev', nlev)
  call vmix_put(Vmix_params,  'prandtl',  0.0d0)
  call vmix_put(Vmix_vars, 'nlev',     nlev)
!  call vmix_put(Vmix_vars, 'visc',     0.0d0)
  Vmix_vars%visc_iface => viscosity
!  call vmix_put(Vmix_vars, 'diff',     0.0d0)
  Vmix_vars%diff_iface => diffusivity
!  call vmix_put(Vmix_vars, 'zw',       iface_depth)
  Vmix_vars%z_iface => iface_depth

  select case (trim(mixing))
    case ('BryanLewis')
      ! Bryan-Lewis mixing
!      print*, "Reading in Bryan-Lewis namelist!"
      read(*, nml=BryanLewis_nml)

      call vmix_init_bkgnd(Vmix_vars, Vmix_params, Vmix_BL_params, &
                           1, 1, vdc1, vdc2, linv, dpth)
      call vmix_coeffs_bkgnd(Vmix_vars, Vmix_BL_params, 1)
    case default
      print*, trim(mixing), " is not a valid mixing choice at this time!"
  end select
  
  ! Get viscosity and diffusivity out of CVMix datatypes
  ! (Not needed if you use pointers instead of vmix_put!)
!  diffusivity = Vmix_vars%diff_iface(:,:)
!  viscosity   = Vmix_vars%visc_iface(:)
  
  do kw=1,nlev+1
    print*, iface_depth(kw), diffusivity(kw,1)
  end do

!EOC

End program vmix_driver
