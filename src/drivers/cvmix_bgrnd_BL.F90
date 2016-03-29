!BOP
!\newpage
! !ROUTINE: cvmix_BL_driver

! !DESCRIPTION: A routine to test the Bryan-Lewis implementation of static
!  background mixing. Inputs are BL coefficients in two columns, one that
!  represents tropical latitudes and one that represents subtropical
!  latitudes. All memory is declared in the driver, and the CVMix data type
!  points to the local variables.
!\\
!\\

! !INTERFACE:

Subroutine cvmix_BL_driver(nlev, max_nlev, ocn_depth)

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                 &
                                    cvmix_zero,               &
                                    cvmix_data_type,          &
                                    cvmix_global_params_type
  use cvmix_background,      only : cvmix_init_bkgnd,         &
                                    cvmix_coeffs_bkgnd,       &
                                    cvmix_get_bkgnd_real_2D,  &
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
  integer, intent(in)        :: nlev,        &! number of levels for column
                                max_nlev      ! number of columns in memory
  real(cvmix_r8), intent(in) :: ocn_depth   ! Depth of ocn

!EOP
!BOC

  ! Global parameter
  integer, parameter :: ncol = 2

  ! CVMix datatypes
  type(cvmix_data_type)         , dimension(ncol) :: CVmix_vars_pointer,      &
                                                     CVmix_vars_memcopy
  type(cvmix_global_params_type)                  :: CVmix_params
  ! Column 1 uses the params saved in module, Column 2 uses this one
  type(cvmix_bkgnd_params_type)                :: CVmix_BL_params

  ! Will use 2 columns, diffusivities be 2 x nlev+1
  ! iface_depth is the depth of each interface;  same in both columns
  real(cvmix_r8), dimension(:),   allocatable, target :: iface_depth
  real(cvmix_r8), dimension(:,:), allocatable, target :: Mdiff, Tdiff

  ! Namelist variables
  ! Bryan-Lewis mixing parameters for column 1
  real(cvmix_r8) :: col1_vdc1, col1_vdc2, col1_linv, col1_dpth
  ! Bryan-Lewis mixing parameters for column 2
  real(cvmix_r8) :: col2_vdc1, col2_vdc2, col2_linv, col2_dpth

  ! array indices
  integer :: i, icol,kw
  ! file indices
#ifdef _NETCDF
  integer, dimension(2) :: fid
#else
  integer, dimension(6) :: fid
#endif

  ! Namelists that may be read in, depending on desired mixing scheme
  namelist/BryanLewis1_nml/col1_vdc1, col1_vdc2, col1_linv, col1_dpth
  namelist/BryanLewis2_nml/col2_vdc1, col2_vdc2, col2_linv, col2_dpth

  print*, "Active levels: ", nlev
  print*, "Levels allocated in memory: ", max_nlev

  ! Calculate depth of cell interfaces based on number of levels and ocean
  ! depth (also allocate memory for diffusivity and viscosity)
  allocate(iface_depth(max_nlev+1))
  iface_depth(1) = cvmix_zero

  ! Depth is 0 at sea level and negative at ocean bottom in CVMix
  do kw = 2,max_nlev+1
    iface_depth(kw) = iface_depth(kw-1) - ocn_depth/real(nlev,cvmix_r8)
  end do

  ! Allocate memory to store viscosity and diffusivity values (for pointer)
  allocate(Mdiff(ncol,max_nlev+1), Tdiff(ncol,max_nlev+1))

  ! Initialization for CVMix data types
  call cvmix_put(CVmix_params,  'max_nlev', max_nlev)
  call cvmix_put(CVmix_params,  'prandtl',  cvmix_zero)
  do icol=1,ncol
    CVmix_vars_pointer(icol)%nlev=nlev
    CVmix_vars_pointer(icol)%max_nlev=max_nlev
    ! Point CVmix_vars_pointer values to memory allocated above
    CVmix_vars_pointer(icol)%Mdiff_iface => Mdiff(icol,:)
    CVmix_vars_pointer(icol)%Tdiff_iface => Tdiff(icol,:)
    CVmix_vars_pointer(icol)%zw_iface    => iface_depth

    ! Copy values into CVmix_vars_memcopy
    call cvmix_put(CVmix_vars_memcopy(icol), 'nlev',     nlev)
    call cvmix_put(CVmix_vars_memcopy(icol), 'max_nlev', max_nlev)
    call cvmix_put(CVmix_vars_memcopy(icol), 'Mdiff',    cvmix_zero)
    call cvmix_put(CVmix_vars_memcopy(icol), 'Tdiff',    cvmix_zero)
    call cvmix_put(CVmix_vars_memcopy(icol), 'zw_iface', iface_depth)
  end do

  ! Read B-L parameters for columns
  read(*, nml=BryanLewis1_nml)
  read(*, nml=BryanLewis2_nml)

  ! Two different columns are tested with each of two different memory options:
  ! Column 1 uses the internal background parameter dataset while column 2 uses
  ! CVmix_BL_params; for the memory copy, column 1 tests the routine
  ! cvmix_get_bkgnd_real_2D() rather than cvmix_coeffs_bkgnd.

  ! Pointer test
  call cvmix_init_bkgnd(CVmix_vars_pointer(1), col1_vdc1, col1_vdc2,          &
                        col1_linv, col1_dpth, CVmix_params)
  call cvmix_init_bkgnd(CVmix_vars_pointer(2), col2_vdc1, col2_vdc2,          &
                        col2_linv, col2_dpth, CVMix_params,                   &
                        CVmix_bkgnd_params_user=CVmix_BL_params)
  call cvmix_coeffs_bkgnd(CVmix_vars_pointer(1))
  call cvmix_coeffs_bkgnd(CVmix_vars_pointer(2),                              &
                          CVmix_bkgnd_params_user=CVmix_BL_params)

  ! Memcopy test
  call cvmix_init_bkgnd(CVmix_vars_memcopy(1), col1_vdc1, col1_vdc2,          &
                        col1_linv, col1_dpth, CVmix_params)
  call cvmix_init_bkgnd(CVmix_vars_memcopy(2), col2_vdc1, col2_vdc2,          &
                        col2_linv, col2_dpth, CVMix_params,                   &
                        CVmix_bkgnd_params_user=CVmix_BL_params)
  CVmix_vars_memcopy(1)%Mdiff_iface = reshape(                                &
                            cvmix_get_bkgnd_real_2D('static_Mdiff'),(/nlev+1/))
  CVmix_vars_memcopy(1)%Tdiff_iface = reshape(                                &
                            cvmix_get_bkgnd_real_2D('static_Tdiff'),(/nlev+1/))
  call cvmix_coeffs_bkgnd(CVmix_vars_memcopy(2),                              &
                          CVmix_bkgnd_params_user=CVmix_BL_params)

  ! Output
#ifdef _NETCDF
  ! data will have diffusivity from both columns (needed for NCL script)
  call cvmix_io_open(fid(1), "data_pointer.nc", "nc")
  call cvmix_io_open(fid(2), "data_memcopy.nc", "nc")
  ! Note: all entries in string of variables to output must be
  !       the same length... hence the space in "Tdiff"
  call cvmix_output_write(fid(1), CVmix_vars_pointer,                         &
                          (/"zw_iface", "Tdiff   "/))
  call cvmix_output_write(fid(2), CVmix_vars_memcopy,                         &
                          (/"zw_iface", "Tdiff   "/))
  do i=1,2
    call cvmix_output_write_att(fid(i), "long_name", "tracer diffusivity",    &
                                var_name="Tdiff")
    call cvmix_output_write_att(fid(i), "units", "m^2/s", var_name="Tdiff")
    call cvmix_output_write_att(fid(i), "long_name", "depth to interface",    &
                                var_name="zw")
    call cvmix_output_write_att(fid(i), "positive", "up", var_name="zw")
    call cvmix_output_write_att(fid(i), "units", "m", var_name="zw")
    call cvmix_io_close(fid(i))
  end do
#else
  ! data will have diffusivity from both columns (needed for NCL script)
  call cvmix_io_open(fid(1), "data_pointer.out", "ascii")
  call cvmix_io_open(fid(2), "data_memcopy.out", "ascii")
  ! col1 will just have diffusivity from low lat
  call cvmix_io_open(fid(3), "col1_pointer.out", "ascii")
  call cvmix_io_open(fid(4), "col1_memcopy.out", "ascii")
  ! col2 will just have diffusivity from high lat
  call cvmix_io_open(fid(5), "col2_pointer.out", "ascii")
  call cvmix_io_open(fid(6), "col2_memcopy.out", "ascii")

  ! Note: all entries in string of variables to output must be
  !       the same length... hence the space in "Tdiff"
  call cvmix_output_write(fid(1), CVmix_vars_pointer,                         &
                          (/"zw_iface", "Tdiff   "/))
  call cvmix_output_write(fid(2), CVmix_vars_memcopy,                         &
                          (/"zw_iface", "Tdiff   "/))
  call cvmix_output_write(fid(3), CVmix_vars_pointer(1),                      &
                          (/"zw_iface", "Tdiff   "/))
  call cvmix_output_write(fid(4), CVmix_vars_memcopy(1),                      &
                          (/"zw_iface", "Tdiff   "/))
  call cvmix_output_write(fid(5), CVmix_vars_pointer(2),                      &
                          (/"zw_iface", "Tdiff   "/))
  call cvmix_output_write(fid(6), CVmix_vars_memcopy(2),                      &
                          (/"zw_iface", "Tdiff   "/))

  do i=1,6
    call cvmix_io_close(fid(i))
  end do
#endif

!EOC

End Subroutine cvmix_BL_driver
