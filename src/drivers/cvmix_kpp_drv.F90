!BOP
!\newpage
! !ROUTINE: cvmix_kpp_driver

! !DESCRIPTION: A routine to test the KPP module.
!\\
!\\

! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !INTERFACE:

Subroutine cvmix_kpp_driver(nlev)

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                 &
                                    cvmix_strlen,             &
                                    cvmix_data_type
  use cvmix_kpp,             only : cvmix_init_kpp,                           &
                                    cvmix_put_kpp,                            &
                                    cvmix_kpp_compute_OBL_depth,              &
                                    cvmix_kpp_compute_turbulent_scales,       &
                                    cvmix_kpp_compute_shape_function_coeffs,  &
                                    cvmix_coeffs_kpp
  use cvmix_put_get,         only : cvmix_put
  use cvmix_io,              only : cvmix_io_open,            &
                                    cvmix_output_write,       &
                                    cvmix_output_write_att,   &
                                    cvmix_io_close

  Implicit None

! !INPUT PARAMETERS:
  integer, intent(in) :: nlev

!EOP
!BOC

  ! CVMix datatypes
  type(cvmix_data_type)       :: CVmix_vars

  real(cvmix_r8), dimension(:,:), allocatable, target :: diffusivity
  real(cvmix_r8), dimension(:),   allocatable, target :: viscosity
  real(cvmix_r8), dimension(:),   allocatable, target :: zt, zw_iface, Ri_bulk
  real(cvmix_r8), dimension(:),   allocatable, target :: w_m, w_s, zeta
  real(cvmix_r8), dimension(:,:), allocatable, target :: TwoDArray
  real(cvmix_r8), dimension(4) :: shape_coeffs
  integer :: fid, kt, kw, nlev2
  real(cvmix_r8) :: hmix, ri_crit, layer_thick
  character(len=cvmix_strlen) :: interp_type

  namelist/kpp_col1_nml/hmix, ri_crit, layer_thick, interp_type

  print*, "Test 1: determining OBL depth"
  print*, "----------"

  hmix = -15.0_cvmix_r8
  ri_crit = 0.3_cvmix_r8
  interp_type = 'quadratic'
  read(*, nml=kpp_col1_nml)

  allocate(diffusivity(nlev+1,2))
  diffusivity = 0.0_cvmix_r8
  diffusivity(2,1) = 10.0_cvmix_r8
  diffusivity(3,1) = 5.0_cvmix_r8
  diffusivity(4,1) = 1.0_cvmix_r8
  allocate(viscosity(nlev+1))
  allocate(zt(nlev), zw_iface(nlev+1), Ri_bulk(nlev))
  do kw=1,nlev+1
    zw_iface(kw) = -layer_thick*real(kw-1, cvmix_r8)
  end do
  do kt=1,nlev
    zt(kt) = 0.5_cvmix_r8*(zw_iface(kt)+zw_iface(kt+1))
    if (zw_iface(kt+1).gt.hmix) then
      Ri_bulk(kt) = 0.0_cvmix_r8
    else
      if (Ri_bulk(kt-1).eq.0) then
        ! Exact integration for average value over first cell with non-zero
        ! Ri_bulk
        Ri_bulk(kt) = 0.25_cvmix_r8*ri_crit*(zw_iface(kt+1)-hmix)**2/layer_thick
      else
        Ri_bulk(kt) = 0.5_cvmix_r8*ri_crit*(hmix-zt(kt))
      end if
    end if
  end do

  call cvmix_put(CVmix_vars, 'nlev', nlev)
  call cvmix_put(CVmix_vars, 'ocn_depth', layer_thick*real(nlev,cvmix_r8))
  call cvmix_put(CVmix_vars, 'surf_fric', 1.0_cvmix_r8)
  call cvmix_put(CVmix_vars, 'surf_buoy', 100.0_cvmix_r8)
  call cvmix_put(CVmix_vars, 'Coriolis', 1e-4_cvmix_r8)
  CVmix_vars%diff_iface => diffusivity(:,:)
  CVmix_vars%visc_iface => viscosity(:)
  CVmix_vars%zt         => zt(:)
  CVmix_vars%zw_iface   => zw_iface(:)
  CVmix_vars%Rib        => Ri_bulk(:)

  call cvmix_init_kpp(ri_crit=ri_crit, vonkarman=0.4_cvmix_r8,                &
                      interp_type=interp_type)
  call cvmix_kpp_compute_OBL_depth(CVmix_vars)
  print*, "OBL depth = ", CVmix_vars%OBL_depth
  print*, "kt of cell containing OBL depth = ", CVmix_vars%kOBL_depth
  call cvmix_coeffs_kpp(CVmix_vars)

#ifdef _NETCDF
  call cvmix_io_open(fid, "data.nc", "nc")
#else
  call cvmix_io_open(fid, "data.out", "ascii")
#endif

  call cvmix_output_write(fid, CVmix_vars, (/"zt     ", "zw     ", "Ri_bulk", &
                                             "diff   "/))
#ifdef _NETCDF
  call cvmix_output_write_att(fid, "Interpolation", interp_type)
  call cvmix_output_write_att(fid, "analytic_OBL_depth", hmix-2.0_cvmix_r8)
  call cvmix_output_write_att(fid, "computed_OBL_depth", CVmix_vars%OBL_depth)
  call cvmix_output_write_att(fid, "kOBL_depth", CVmix_vars%kOBL_depth)
#endif

  call cvmix_io_close(fid)

  print*, ""
  print*, "Test 2: Computing G(sigma)"
  print*, "----------"
  call cvmix_kpp_compute_shape_function_coeffs(0.0_cvmix_r8, 0.0_cvmix_r8,    &
                                               shape_coeffs)
  print*, "Coefficients are: "
  print*, shape_coeffs(1), shape_coeffs(2), shape_coeffs(3), shape_coeffs(4) 

  print*, ""
  print*, "Test 3: determining phi_m and phi_s (inversely proportional to ",  &
          "w_m and w_s, respectively)"
  print*, "----------"
  call cvmix_put_kpp('vonkarman', 1.0_cvmix_r8)
  nlev2 = 220
  allocate(w_m(nlev2+1), w_s(nlev2+1), zeta(nlev2+1))
  ! Note: zeta = sigma*OBL_depth/MoninObukhov constant
  !       zeta < 0 => unstable flow
  !       zeta > 0 => stable flow
  do kw=1, nlev2+1
    zeta(kw) = -2.0_cvmix_r8 + 2.2_cvmix_r8*real(kw-1,cvmix_r8)/real(nlev2,cvmix_r8)
  end do
  ! Typically the first argument of compute_turbulent_scales is sigma, and then
  ! the routine calculates zeta based on the next three parameters. Setting
  ! OBL_depth = surf_buoy_force = surf_fric_vel = 1 (with von Karman = 1 as
  ! well) => sigma = zeta
  call cvmix_kpp_compute_turbulent_scales(zeta, 1.0_cvmix_r8, 1.0_cvmix_r8,  &
                                          1.0_cvmix_r8, w_m, w_s)

  allocate(TwoDArray(nlev2+1,3))
  TwoDArray(:,1) = zeta
  TwoDArray(:,2) = 1.0_cvmix_r8/w_m ! phi_m
  TwoDArray(:,3) = 1.0_cvmix_r8/w_s ! phi_s
#ifdef _NETCDF
  call cvmix_io_open(fid, "test3.nc", "nc")
#else
  call cvmix_io_open(fid, "test3.out", "ascii")
#endif
  call cvmix_output_write(fid, "data", (/"nrow", "ncol"/), TwoDArray)
  call cvmix_io_close(fid)
#ifdef _NETCDF
  print*, "Done! Data is stored in test3.nc, run plot_flux_profiles.ncl to see output."
#else
  print*, "Done! Data is stored in test3.out, run plot_flux_profiles.ncl to see output."
#endif
  print*, ""
  deallocate(TwoDArray)
  deallocate(zeta, w_m, w_s)

!EOC

End Subroutine cvmix_kpp_driver
