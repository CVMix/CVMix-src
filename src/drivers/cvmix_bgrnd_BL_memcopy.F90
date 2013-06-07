!BOP
!\newpage
! !ROUTINE: cvmix_BL_memcopy_driver

! !DESCRIPTION: A routine to test the Bryan-Lewis implementation of static
!  background mixing. Inputs are BL coefficients in two columns, one that
!  represents tropical latitudes and one that represents subtropical
!  latitudes. All memory is declared in the driver and then copied into the
!  CVMix data structures.
!\\
!\\

! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !INTERFACE:

Subroutine cvmix_BL_memcopy_driver(nlev, ocn_depth)

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                 &
                                    cvmix_data_type,          &
                                    cvmix_global_params_type
  use cvmix_background,      only : cvmix_init_bkgnd,         &
                                    cvmix_coeffs_bkgnd,       &
                                    cvmix_bkgnd_params_type
  use cvmix_put_get,         only : cvmix_put
  use cvmix_io,              only : cvmix_io_open,            &
                                    cvmix_output_write,       &
#ifdef _NETCDF
                                    cvmix_output_write_att,   &
#endif
                                    cvmix_io_close

  Implicit None
  
! !INPUT PARAMETERS:
  integer, intent(in)        :: nlev        ! number of levels for column
  real(cvmix_r8), intent(in) :: ocn_depth   ! Depth of ocn

!EOP
!BOC

! !LOCAL VARIABLES:

  ! Global parameter
  integer, parameter :: ncol = 2

  ! CVMix datatypes
  type(cvmix_data_type)         , dimension(2) :: CVmix_vars
  type(cvmix_global_params_type)               :: CVmix_params
  type(cvmix_bkgnd_params_type) , dimension(2) :: CVmix_BL_params

  ! iface_depth is the depth of each interface;  same in both columns
  real(cvmix_r8), dimension(:), allocatable :: iface_depth

  ! Namelist variables
  ! Bryan-Lewis mixing parameters for column 1
  real(cvmix_r8) :: col1_vdc1, col1_vdc2, col1_linv, col1_dpth
  ! Bryan-Lewis mixing parameters for column 2
  real(cvmix_r8) :: col2_vdc1, col2_vdc2, col2_linv, col2_dpth

  ! array indices
  integer :: icol,kw
  ! file indices
#ifdef _NETCDF
  integer :: fid
#else
  integer :: fid1, fid2, fid3
#endif

  ! Namelists that may be read in, depending on desired mixing scheme
  namelist/BryanLewis1_nml/col1_vdc1, col1_vdc2, col1_linv, col1_dpth
  namelist/BryanLewis2_nml/col2_vdc1, col2_vdc2, col2_linv, col2_dpth

  ! Calculate depth of cell interfaces based on number of levels and ocean
  ! depth
  allocate(iface_depth(nlev+1))
  iface_depth(1) = 0.0_cvmix_r8
  
  ! Depth is 0 at sea level and negative at ocean bottom in CVMix
  do kw = 2,nlev+1
    iface_depth(kw) = iface_depth(kw-1) - ocn_depth/real(nlev,cvmix_r8)
  end do

  ! Initialization for CVMix data types
  ! Note that cvmix_put will allocate memory for arrays as needed
  call cvmix_put(CVmix_params, 'max_nlev', nlev)
  call cvmix_put(CVmix_params, 'prandtl',  0.0_cvmix_r8)
  do icol=1,2
    call cvmix_put(CVmix_vars(icol), 'nlev',     nlev)
    call cvmix_put(CVmix_vars(icol), 'visc',     0.0_cvmix_r8)
    call cvmix_put(CVmix_vars(icol), 'diff',     0.0_cvmix_r8)
    call cvmix_put(CVmix_vars(icol), 'zw_iface', iface_depth)
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
  call cvmix_io_open(fid, "data.nc", "nc")
  ! Note: all entries in string of variables to output must be
  !       the same length... hence the space in "diff "
  call cvmix_output_write(fid, CVmix_vars, (/"depth", "diff "/))
  call cvmix_output_write_att(fid, "long_name", "tracer diffusivity",         &
                              var_name="diff")
  call cvmix_output_write_att(fid, "units", "m^2/s", var_name="diff")
  call cvmix_output_write_att(fid, "long_name", "depth to interface",         &
                              var_name="depth")
  call cvmix_output_write_att(fid, "positive", "up", var_name="depth")
  call cvmix_output_write_att(fid, "units", "m", var_name="depth")
  call cvmix_io_close(fid)
#else
  ! data will have diffusivity from both columns (needed for NCL script)
  call cvmix_io_open(fid1, "data.out", "ascii")
  ! col1 will just have diffusivity from low lat
  call cvmix_io_open(fid2, "col1.out", "ascii")
  ! col2 will just have diffusivity from high lat
  call cvmix_io_open(fid3, "col2.out", "ascii")

  ! Note: all entries in string of variables to output must be
  !       the same length... hence the space in "diff "
  call cvmix_output_write(fid1, CVmix_vars,    (/"depth", "diff "/))
  call cvmix_output_write(fid2, CVmix_vars(1), (/"depth", "diff "/))
  call cvmix_output_write(fid3, CVmix_vars(2), (/"depth", "diff "/))

  call cvmix_io_close(fid1)
  call cvmix_io_close(fid2)
  call cvmix_io_close(fid3)
#endif

!EOC

End Subroutine cvmix_BL_memcopy_driver
