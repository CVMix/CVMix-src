!BOP

! !ROUTINE: cvmix_BL_driver_pointers

! !DESCRIPTION: A stand-alone driver for the CVMix package. This particular
!  driver generates the Bryan-Lewis coefficients in two columns, one that
!  represents tropical latitudes and one that represents subtropical
!  latitudes. All memory is declared in the driver, and the CVMix data type
!  points to the local variables.
!\\
!\\

! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !INTERFACE:

Program cvmix_BL_driver_pointers

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                 &
                                    cvmix_data_type,          &
                                    cvmix_global_params_type, &
                                    cvmix_bkgnd_params_type
  use cvmix_background,      only : cvmix_init_bkgnd,         &
                                    cvmix_coeffs_bkgnd
  use cvmix_put_get,         only : cvmix_put
  use cvmix_output,          only : cvmix_output_open,        &
                                    cvmix_output_write,       &
                                    cvmix_output_close
!EOP
!BOC

  Implicit None

  type(cvmix_data_type)         , dimension(2) :: CVmix_vars
  type(cvmix_global_params_type)               :: CVmix_params
  type(cvmix_bkgnd_params_type) , dimension(2) :: CVmix_BL_params

  ! Will use 2 columns, viscosity will be 2 x nlev+1 and diffusivity will 
  ! be 2 x nlev+1 x 1 (diffusivity is 2D in column)
  ! iface_depth is the depth of each interface;  same in both columns
  real(cvmix_r8), dimension(:),     allocatable, target :: iface_depth
  real(cvmix_r8), dimension(:,:),   allocatable, target :: viscosity
  real(cvmix_r8), dimension(:,:,:), allocatable, target :: diffusivity

  ! Global parameter
  integer, parameter :: ncol = 2
  ! array indices
  integer :: icol,kw
  ! file indices
#ifdef _NETCDF
  integer :: fid
#else
  integer :: fid1, fid2, fid3
#endif

  ! Namelist variables
  ! 1) General mixing parameters
  integer                :: nlev      ! number of levels for column
  real(cvmix_r8)         :: ocn_depth ! Depth of ocn
  ! 2) Bryan-Lewis mixing parameters for column 1
  real(cvmix_r8)         :: col1_vdc1, col1_vdc2, col1_linv, col1_dpth
  ! 3) Bryan-Lewis mixing parameters for column 2
  real(cvmix_r8)         :: col2_vdc1, col2_vdc2, col2_linv, col2_dpth

  ! Namelists that may be read in, depending on desired mixing scheme
  namelist/cvmix_nml/nlev, ocn_depth
  namelist/BryanLewis1_nml/col1_vdc1, col1_vdc2, col1_linv, col1_dpth
  namelist/BryanLewis2_nml/col2_vdc1, col2_vdc2, col2_linv, col2_dpth

  ! Read general mixing parameters
  read(*, nml=cvmix_nml)

  ! Calculate depth of cell interfaces based on number of levels and ocean
  ! depth (also allocate memory for diffusivity and viscosity)
  allocate(iface_depth(nlev+1))
  iface_depth(1) = 0.0_cvmix_r8
  do kw = 2,nlev+1
    iface_depth(kw) = iface_depth(kw-1) + ocn_depth/dble(nlev)
  end do

  ! Allocate memory to store viscosity and diffusivity values
  allocate(diffusivity(2,nlev+1,1), viscosity(2,nlev+1))

  ! Initialization for CVMix data types
  call cvmix_put(CVmix_params,  'max_nlev', nlev)
  call cvmix_put(CVmix_params,  'prandtl',  0.0_cvmix_r8)
  do icol=1,2
    call cvmix_put(CVmix_vars(icol), 'nlev',     nlev)
    ! Point CVmix_vars values to memory allocated above
    CVmix_vars(icol)%visc_iface => viscosity(icol,:)
    CVmix_vars(icol)%diff_iface => diffusivity(icol,:,:)
    CVmix_vars(icol)%z_iface => iface_depth
  end do

  ! Read / set B-L parameters for column 1
  read(*, nml=BryanLewis1_nml)
  call cvmix_init_bkgnd(CVmix_vars(1), CVmix_params, CVmix_BL_params(1), &
                        col1_vdc1, col1_vdc2, col1_linv, col1_dpth)
  call cvmix_coeffs_bkgnd(CVmix_vars(1), CVmix_BL_params(1), 1)
  
  ! Read / set B-L parameters for column 2
  read(*, nml=BryanLewis2_nml)
  call cvmix_init_bkgnd(CVmix_vars(2), CVmix_params, CVmix_BL_params(2), &
                        col2_vdc1, col2_vdc2, col2_linv, col2_dpth)
  call cvmix_coeffs_bkgnd(CVmix_vars(2), CVmix_BL_params(2), 1)
  
  ! Output
#ifdef _NETCDF
  ! data will have diffusivity from both columns (needed for NCL script)
  call cvmix_output_open(fid, "data.nc", "nc")
  ! Note: all entries in string of variables to output must be
  !       the same length... hence the space in "diff "
  call cvmix_output_write(fid, CVmix_vars, (/"depth", "diff "/))
  call cvmix_output_close(fid)
#else
  ! data will have diffusivity from both columns (needed for NCL script)
  call cvmix_output_open(fid1, "data.out", "ascii")
  ! col1 will just have diffusivity from low lat
  call cvmix_output_open(fid2, "col1.out", "ascii")
  ! col2 will just have diffusivity from high lat
  call cvmix_output_open(fid3, "col2.out", "ascii")

  ! Note: all entries in string of variables to output must be
  !       the same length... hence the space in "diff "
  call cvmix_output_write(fid1, CVmix_vars,    (/"depth", "diff "/))
  call cvmix_output_write(fid2, CVmix_vars(1), (/"depth", "diff "/))
  call cvmix_output_write(fid3, CVmix_vars(2), (/"depth", "diff "/))

  call cvmix_output_close(fid1)
  call cvmix_output_close(fid2)
  call cvmix_output_close(fid3)
#endif

!EOC

End program cvmix_BL_driver_pointers
