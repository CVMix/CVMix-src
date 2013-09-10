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

Subroutine cvmix_kpp_driver()

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                 &
                                    cvmix_strlen,             &
                                    cvmix_data_type
  use cvmix_kpp,             only : cvmix_init_kpp,                           &
                                    cvmix_put_kpp,                            &
                                    cvmix_kpp_compute_OBL_depth,              &
                                    cvmix_kpp_compute_kOBL_depth,             &
                                    cvmix_kpp_compute_bulk_Richardson,        &
                                    cvmix_kpp_compute_unresolved_shear,       &
                                    cvmix_kpp_compute_turbulent_scales,       &
                                    cvmix_kpp_compute_shape_function_coeffs,  &
                                    cvmix_coeffs_kpp
  use cvmix_put_get,         only : cvmix_put
  use cvmix_io,              only : cvmix_io_open,            &
                                    cvmix_output_write,       &
                                    cvmix_output_write_att,   &
                                    cvmix_io_close

  Implicit None

!EOP
!BOC

  ! CVMix datatypes
  type(cvmix_data_type)       :: CVmix_vars1, CVmix_vars4

  real(cvmix_r8), dimension(:,:), allocatable, target :: diffusivity
  real(cvmix_r8), dimension(:),   allocatable, target :: viscosity
  real(cvmix_r8), dimension(:),   allocatable, target :: zt, zw_iface, Ri_bulk
  real(cvmix_r8), dimension(:),   allocatable, target :: w_m, w_s, zeta
  real(cvmix_r8), dimension(:,:), allocatable, target :: TwoDArray
  real(cvmix_r8), dimension(:),   allocatable, target :: buoyancy, shear_sqr, &
                                                         delta_vel_sqr,       &
                                                         buoy_freq
  real(cvmix_r8), dimension(:,:), allocatable, target :: hor_vel
  real(cvmix_r8), dimension(2)                        :: ref_vel
  real(cvmix_r8), dimension(4) :: shape_coeffs
  integer :: fid, kt, kw, nlev1, nlev3, nlev4, nlev5
  real(cvmix_r8) :: hmix1, hmix5, ri_crit, layer_thick, OBL_depth, N, Nsqr
  real(cvmix_r8) :: kOBL_depth, Bslope, Vslope
  character(len=cvmix_strlen) :: interp_type_t1, interp_type_t4
  logical :: ltest1, ltest2, ltest3, ltest4, ltest5 ! True => run specfied test
  logical :: lnoDGat1 ! True => G'(1) = 0 (in test 4)

  namelist/kpp_col1_nml/ltest1, nlev1, layer_thick, interp_type_t1, hmix1,    &
                        ri_crit
  namelist/kpp_col2_nml/ltest2
  namelist/kpp_col3_nml/ltest3, nlev3
  namelist/kpp_col4_nml/ltest4, interp_type_t4, OBL_depth, lnoDGat1

  ! Read namelists

  ! Defaults for test 1 
  ltest1         = .false.
  nlev1          = 4
  layer_thick    = 10.0_cvmix_r8
  hmix1          = -15.0_cvmix_r8
  ri_crit        = 0.3_cvmix_r8
  interp_type_t1 = 'quadratic'

  ! Defaults for test 2 
  ltest2 = .false.

  ! Defaults for test 3 
  ltest3 = .false.
  nlev3  = 220

  ! Defaults for test 4 
  ltest4         = .false.
  OBL_depth      = 14.0_cvmix_r8 
  interp_type_t4 = 'quadratic'
  lnoDGat1       = .true.

  ! Defaults for test 5
  ltest5 = .true.

  read(*, nml=kpp_col1_nml)
  read(*, nml=kpp_col2_nml)
  read(*, nml=kpp_col3_nml)
  read(*, nml=kpp_col4_nml)

  ! Test 1: user sets up levels via namelist (constant thickness) and specifies
  !         critical Richardson number as well as depth parameter hmix1. The
  !         bulk Richardson number is assumed to be 0 from surface to hmix1 and
  !         then increases linearly at a rate of Ri_crit/2 (so bulk Richardson
  !         number = Ri_crit at hmix1+2). For computation, though, the average
  !         bulk Richardson number (integral over vertical layer divided by
  !         layer thickness) is stored at cell centers and then interpolated
  !         (user can specify linear, quadratic or cubic interpolant) between
  !         cell centers. OBL_depth is set to depth where interpolated bulk
  !         Richardson number = Ri_crit; level-center depth (zt) and bulk 
  !         Richardson numbers are written out to test1.nc or test1.out
  if (ltest1) then
    print*, "Test 1: determining OBL depth"
    print*, "----------"

    ! Initialize parameter datatype and set up column
    call cvmix_init_kpp(ri_crit=ri_crit, interp_type=interp_type_t1)
    call cvmix_put(CVmix_vars1, 'nlev', nlev1)
    call cvmix_put(CVmix_vars1, 'ocn_depth', layer_thick*real(nlev1,cvmix_r8))

    ! Set up vertical levels (centers and interfaces) and compute bulk
    ! Richardson number
    allocate(zt(nlev1), zw_iface(nlev1+1), Ri_bulk(nlev1))
    do kw=1,nlev1+1
      zw_iface(kw) = -layer_thick*real(kw-1, cvmix_r8)
    end do
    do kt=1,nlev1
      zt(kt) = 0.5_cvmix_r8*(zw_iface(kt)+zw_iface(kt+1))
      if (zw_iface(kt+1).gt.hmix1) then
        Ri_bulk(kt) = 0.0_cvmix_r8
      else
        if (Ri_bulk(kt-1).eq.0) then
          ! Exact integration for average value over first cell with non-zero
          ! Ri_bulk
          Ri_bulk(kt) = 0.25_cvmix_r8*ri_crit*(zw_iface(kt+1)-hmix1)**2 /     &
                        layer_thick
        else
          Ri_bulk(kt) = 0.5_cvmix_r8*ri_crit*(hmix1-zt(kt))
        end if
      end if
    end do

    CVmix_vars1%zt       => zt(:)
    CVmix_vars1%zw_iface => zw_iface(:)
    CVmix_vars1%Rib      => Ri_bulk(:)

    ! Compute OBL depth
    call cvmix_kpp_compute_OBL_depth(CVmix_vars1)

    ! Output to screen and file
    print*, "OBL depth = ", CVmix_vars1%OBL_depth
    print*, "kw of interface above OBL depth = ", floor(CVmix_vars1%kOBL_depth)
    print*, "kt of cell center above OBL depth = ", nint(CVmix_vars1%kOBL_depth)-1

#ifdef _NETCDF
    call cvmix_io_open(fid, "test1.nc", "nc")
#else
    call cvmix_io_open(fid, "test1.out", "ascii")
#endif

    call cvmix_output_write(fid, CVmix_vars1, (/"zt     ", "zw     ",         &
                                                "Ri_bulk"/))
#ifdef _NETCDF
    call cvmix_output_write_att(fid, "Interpolation", interp_type_t1)
    call cvmix_output_write_att(fid, "analytic_OBL_depth", -hmix1+2.0_cvmix_r8)
    call cvmix_output_write_att(fid, "computed_OBL_depth",                    &
                                CVmix_vars1%OBL_depth)
    call cvmix_output_write_att(fid, "kOBL_depth", CVmix_vars1%kOBL_depth)
#endif

    call cvmix_io_close(fid)

    deallocate(zt, zw_iface, Ri_bulk)
  end if ! ltest for Test 1

  ! Test 2: Compute coefficients of shape function G(sigma) when G(1) = 0 and
  !         G'(1) = 0. Result should be G(sigma) = sigma - 2sigma^2 + sigma^3
  if (ltest2) then
    print*, ""
    print*, "Test 2: Computing G(sigma)"
    print*, "----------"

    call cvmix_kpp_compute_shape_function_coeffs(0.0_cvmix_r8, 0.0_cvmix_r8,  &
                                                 shape_coeffs)
    write(*,"(1X,A,4F7.3)") "Coefficients are: ", shape_coeffs
  end if ! ltest for test 2

  ! Test 3: Recreate Figure B1 in LMD94 (phi(zeta)). Note that von Karman,
  !         surface buoyancy forcing, and surface velocity are set such that
  !         Monin-Obukhov constant = 1 => zeta = sigma.
  if (ltest3) then
    print*, ""
    print*, "Test 3: determining phi_m and phi_s (inversely proportional to ",&
            "w_m and w_s, respectively)"
    print*, "----------"
    call cvmix_put_kpp('vonkarman', 1.0_cvmix_r8)
    call cvmix_put_kpp('surf_layer_ext', 1.0_cvmix_r8)
    allocate(w_m(nlev3+1), w_s(nlev3+1), zeta(nlev3+1))
    ! Note: zeta = sigma*OBL_depth/MoninObukhov constant
    !       zeta < 0 => unstable flow
    !       zeta > 0 => stable flow
    do kw=1, nlev3+1
      zeta(kw) = -2.0_cvmix_r8 + 2.2_cvmix_r8 * real(kw-1,cvmix_r8) /         &
                 real(nlev3,cvmix_r8)
    end do
    ! Typically the first argument of compute_turbulent_scales is sigma, and then
    ! the routine calculates zeta based on the next three parameters. Setting
    ! OBL_depth = surf_buoy_force = surf_fric_vel = 1 (with von Karman = 1 as
    ! well) => sigma = zeta
    call cvmix_kpp_compute_turbulent_scales(zeta, 1.0_cvmix_r8, 1.0_cvmix_r8, &
                                            1.0_cvmix_r8, w_m, w_s)

    allocate(TwoDArray(nlev3+1,3))
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
    print*, "Done! Data is stored in test3.nc, run plot_flux_profiles.ncl ",  &
            "to see output."
#else
    print*, "Done! Data is stored in test3.out, run plot_flux_profiles.ncl ", &
            "to see output."
#endif
    print*, ""
    deallocate(TwoDArray)
    deallocate(zeta, w_m, w_s)
  endif ! ltest3

  if (ltest4) then
    print*, ""
    print*, "Test 4: Computing Diffusivity in boundary layer"
    print*, "----------"

    nlev4 = 5
    layer_thick = 5

    ! Set up vertical levels (centers and interfaces) and compute bulk
    ! Richardson number
    allocate(zt(nlev4), zw_iface(nlev4+1))
    do kw=1,nlev4+1
      zw_iface(kw) = -layer_thick*real(kw-1, cvmix_r8)
    end do
    do kt=1,nlev4
      zt(kt) = 0.5_cvmix_r8*(zw_iface(kt)+zw_iface(kt+1))
    end do
    CVmix_vars4%zt       => zt(:)
    CVmix_vars4%zw_iface => zw_iface(:)

    ! Set up diffusivities
    allocate(diffusivity(nlev4+1,2), viscosity(nlev4+1))
    diffusivity = 0.0_cvmix_r8
    diffusivity(2,1) = 10.0_cvmix_r8
    diffusivity(3,1) = 5.0_cvmix_r8
    diffusivity(4,1) = 1.0_cvmix_r8
    CVmix_vars4%diff_iface => diffusivity
    CVmix_vars4%visc_iface => viscosity

    ! Set physical properties of column for test 4
    call cvmix_put(CVmix_vars4, 'nlev', nlev4)
    call cvmix_put(CVmix_vars4, 'ocn_depth', layer_thick*real(nlev4,cvmix_r8))
    call cvmix_put(CVmix_vars4, 'OBL_depth', OBL_depth)
    call cvmix_put(CVmix_vars4, 'kOBL_depth',                                 &
         cvmix_kpp_compute_kOBL_depth(zw_iface, zt, OBL_depth))

    print*, "OBL_depth = ", CVmix_vars4%OBL_depth
    print*, "ocean depth = ", CVmix_vars4%ocn_depth
    print*, "kOBL_depth = ", CVmix_vars4%kOBL_depth

    call cvmix_put(CVmix_vars4, 'surf_fric', 1.0_cvmix_r8)
    call cvmix_put(CVmix_vars4, 'surf_buoy', 100.0_cvmix_r8)
    call cvmix_put(CVmix_vars4, 'Coriolis', 1e-4_cvmix_r8)

    call cvmix_init_kpp(ri_crit=ri_crit, vonkarman=0.4_cvmix_r8,              &
                        interp_type2=interp_type_t4, lnoDGat1=lnoDGat1)
    call cvmix_coeffs_kpp(CVmix_vars4)

#ifdef _NETCDF
    call cvmix_io_open(fid, "test4.nc", "nc")
#else
    call cvmix_io_open(fid, "test4.out", "ascii")
#endif

    call cvmix_output_write(fid, CVmix_vars4, (/"zt  ", "zw  ", "diff"/)) 
#ifdef _NETCDF
    call cvmix_output_write_att(fid, "interp_type2", interp_type_t4)
    call cvmix_output_write_att(fid, "OBL_depth", CVmix_vars4%OBL_depth)
#endif

    call cvmix_io_close(fid)

    deallocate(zt, zw_iface)
    deallocate(diffusivity, viscosity)
  end if ! ltest4

  ! Test 5: Recreate figure C1 from LMD94
  if (ltest5) then
    print*, ""
    print*, "Test 5: Computing Bulk Richardson number"
    print*, "----------"

    nlev5 = 10
    layer_thick = 5
    hmix5 = 17
    ! using linear interpolation to match LMD result
    call cvmix_init_kpp(vonkarman=0.4_cvmix_r8, interp_type='linear')

    ! Set up vertical levels (centers and interfaces) and compute bulk
    ! Richardson number
    allocate(zt(nlev5), zw_iface(nlev5+1))
    do kw=1,nlev5+1
      zw_iface(kw) = -layer_thick*real(kw-1, cvmix_r8)
    end do
    do kt=1,nlev5
      zt(kt) = 0.5_cvmix_r8*(zw_iface(kt)+zw_iface(kt+1))
    end do

    ! Compute Br-B(d), |Vr-V(d)|^2, and Vt^2
    allocate(buoyancy(nlev5), delta_vel_sqr(nlev5), hor_vel(nlev5,2),         &
             shear_sqr(nlev5), w_s(nlev5), Ri_bulk(nlev5), buoy_freq(nlev5))

    ref_vel(1) = 0.1_cvmix_r8
    ref_vel(2) = 0.0_cvmix_r8
    N            = 0.01_cvmix_r8
    buoy_freq(:) = N
    Nsqr         = N*N
    Bslope       = -Nsqr
    Vslope       = -0.1_cvmix_r8/real(nlev5*layer_thick-hmix5,cvmix_r8)
    do kt=1,nlev5
      if ((zt(kt).ge.-hmix5).or.(kt.eq.1)) then
        buoyancy(kt)  = Nsqr
        hor_vel(kt,1) = 0.1_cvmix_r8
      else
        if (zw_iface(kt).ge.-hmix5) then
          ! derivatives of buoyancy and horizontal velocity component are
          ! discontinuous in this layer (no change -> non-zero linear change)
          ! so we compute area-average of analytic function over layer
          buoyancy(kt)  = Bslope*(-zw_iface(kt+1)-real(hmix5,cvmix_r8))**2 /  &
                          real(2*layer_thick,cvmix_r8) + Nsqr
          hor_vel(kt,1) = Vslope*(-zw_iface(kt+1)-real(hmix5,cvmix_r8))**2 /  &
                          real(2*layer_thick,cvmix_r8) + 0.1_cvmix_r8
        else
          buoyancy(kt)  = Nsqr+Bslope*(-zt(kt)-real(hmix5,cvmix_r8))
          hor_vel(kt,1) = 0.1_cvmix_r8+Vslope*(-zt(kt)-real(hmix5,cvmix_r8))
        end if
      end if
      ! Compute w_s with zeta=0 per LMD page 393
      ! => w_s = von Karman * surf_fric_vel = 0.4*0.01 = 4e-3
      call cvmix_kpp_compute_turbulent_scales(0.0_cvmix_r8, -zt(kt),          &
                 buoyancy(1), 0.01_cvmix_r8, w_s=w_s(kt))
      hor_vel(kt,2) = 0.0_cvmix_r8
      delta_vel_sqr(kt) = (ref_vel(1)-hor_vel(kt,1))**2 +                     &
                          (ref_vel(2)-hor_vel(kt,2))**2
    end do
!   MNL: tested both interfaces of compute_bulk_Richardson
!    shear_sqr = cvmix_kpp_compute_unresolved_shear(zt, buoy_freq, w_s)
!    Ri_bulk   = cvmix_kpp_compute_bulk_Richardson(zt, (buoyancy(1)-buoyancy), &
!                                                  delta_vel_sqr, shear_sqr)
    Ri_bulk = cvmix_kpp_compute_bulk_Richardson(zt, (buoyancy(1)-buoyancy),   &
                                                delta_vel_sqr, buoy_freq, w_s)
    call cvmix_kpp_compute_OBL_depth(Ri_bulk, zw_iface, OBL_depth, kOBL_depth,&
                                     zt)
    do kt=1,nlev5
      print*, zt(kt), Ri_bulk(kt)
    end do
    print*, "OBL has depth of ", OBL_depth
    deallocate(zt, zw_iface)
    deallocate(buoyancy, delta_vel_sqr, hor_vel, shear_sqr, w_s, Ri_bulk,     &
               buoy_freq)
  end if

!EOC

End Subroutine cvmix_kpp_driver
