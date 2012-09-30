!BOI

! !TITLE: In-code documentation for CVMix
! !AUTHORS: Many contributors from GFDL, LANL, and NCAR
! !AFFILIATION: GFDL, LANL, and NCAR
! !DATE: \today

!EOI
!BOP

! !ROUTINE: vmix_BL_driver_pointers

! !DESCRIPTION: The stand-alone driver for the CVMix package. The type of
!  mixing comes from the cvmix\_nml namelist, and then other namelists
!  contain parameters for the specific mixing methods.
!\\
!\\

! !INTERFACE:

Program vmix_BL_driver_pointers

! !USES:

  use vmix_kinds_and_types
  use vmix_background
  use vmix_convection
  use vmix_put_get
  use vmix_output

!EOP
!BOC

  Implicit None

  type (vmix_data_type)         , dimension(2) :: Vmix_vars
  type (vmix_global_params_type)               :: Vmix_params
  type (vmix_bkgnd_params_type) , dimension(2) :: Vmix_BL_params

  ! Will use 2 columns, viscosity will be 2 x nlev+1 and diffusivity will 
  ! be 2 x nlev+1 x 1 (diffusivity is 2D in column)
  ! iface_depth is the depth of each interface;  same in both columns
  real(kind=vmix_r8), dimension(:),     allocatable, target :: iface_depth
  real(kind=vmix_r8), dimension(:,:),   allocatable, target :: viscosity
  real(kind=vmix_r8), dimension(:,:,:), allocatable, target :: diffusivity

  ! Global parameter
  integer, parameter :: ncol = 2
  ! array indices
  integer :: icol,kw
  ! file indices
  integer :: fid1, fid2, fid3

  ! Namelist variables
  ! 1) General mixing parameters
  integer                    :: nlev      ! number of levels for column
  real(kind=vmix_r8)         :: ocn_depth ! Depth of ocn
  ! 2) Bryan-Lewis mixing parameters for column 1
  real(kind=vmix_r8)         :: col1_vdc1, col1_vdc2, col1_linv, col1_dpth
  ! 3) Bryan-Lewis mixing parameters for column 2
  real(kind=vmix_r8)         :: col2_vdc1, col2_vdc2, col2_linv, col2_dpth

  ! Namelists that may be read in, depending on desired mixing scheme
  namelist/cvmix_nml/nlev, ocn_depth
  namelist/BryanLewis1_nml/col1_vdc1, col1_vdc2, col1_linv, col1_dpth
  namelist/BryanLewis2_nml/col2_vdc1, col2_vdc2, col2_linv, col2_dpth

  ! Read general mixing parameters
  read(*, nml=cvmix_nml)

  ! Calculate depth of cell interfaces based on number of levels and ocean
  ! depth (also allocate memory for diffusivity and viscosity)
  allocate(iface_depth(nlev+1))
  iface_depth(1) = 0.0d0
  do kw = 2,nlev+1
    iface_depth(kw) = iface_depth(kw-1) + ocn_depth/dble(nlev)
  end do

  ! Allocate memory to store viscosity and diffusivity values
  allocate(diffusivity(2,nlev+1,1), viscosity(2,nlev+1))

  ! Initialization for CVMix data types
  call vmix_put(Vmix_params,  'max_nlev', nlev)
  call vmix_put(Vmix_params,  'prandtl',  0.0d0)
  do icol=1,2
    call vmix_put(Vmix_vars(icol), 'nlev',     nlev)
    ! Point Vmix_vars values to memory allocated above
    Vmix_vars(icol)%visc_iface => viscosity(icol,:)
    Vmix_vars(icol)%diff_iface => diffusivity(icol,:,:)
    Vmix_vars(icol)%z_iface => iface_depth
  end do

  ! Read / set B-L parameters for column 1
  read(*, nml=BryanLewis1_nml)
  call vmix_init_bkgnd(Vmix_vars(1), Vmix_params, Vmix_BL_params(1),     &
                       1, 1, col1_vdc1, col1_vdc2, col1_linv, col1_dpth)
  call vmix_coeffs_bkgnd(Vmix_vars(1), Vmix_BL_params(1), 1)
  
  ! Read / set B-L parameters for column 2
  read(*, nml=BryanLewis2_nml)
  call vmix_init_bkgnd(Vmix_vars(2), Vmix_params, Vmix_BL_params(2),     &
                       1, 1, col2_vdc1, col2_vdc2, col2_linv, col2_dpth)
  call vmix_coeffs_bkgnd(Vmix_vars(2), Vmix_BL_params(2), 1)
  
  ! Output (just text for now, will be netCDF soon)
  ! data will have diffusivity from both columns (needed for NCL script)
#ifdef _NETCDF
  call vmix_output_open(fid1, "data.nc", "nc")
#else
  call vmix_output_open(fid1, "data.out", "ascii")
#endif

  ! col1 will just have diffusivity from low lat
#ifdef _NETCDF
  call vmix_output_open(fid2, "col1.nc", "nc")
#else
  call vmix_output_open(fid2, "col1.out", "ascii")
#endif

  ! col2 will just have diffusivity from high lat
#ifdef _NETCDF
  call vmix_output_open(fid3, "col2.nc", "nc")
#else
  call vmix_output_open(fid3, "col2.out", "ascii")
#endif

  call vmix_output_write_diffusivity(fid1, Vmix_vars)
  call vmix_output_write_diffusivity(fid2, Vmix_vars(1))
  call vmix_output_write_diffusivity(fid3, Vmix_vars(2))

  call vmix_output_close(fid1)
  call vmix_output_close(fid2)
  call vmix_output_close(fid3)

!EOC

End program vmix_BL_driver_pointers
