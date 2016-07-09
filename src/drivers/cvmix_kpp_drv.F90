!BOP
!\newpage
! !ROUTINE: cvmix_kpp_driver

! !DESCRIPTION: A routine to test the KPP module.
!\\
!\\

! !INTERFACE:

Subroutine cvmix_kpp_driver()

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                 &
                                    cvmix_zero,               &
                                    cvmix_one,                &
                                    cvmix_strlen,             &
                                    cvmix_data_type
  use cvmix_kpp,             only : cvmix_init_kpp,                           &
                                    cvmix_get_kpp_real,                       &
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
#ifdef _NETCDF
                                    cvmix_output_write_att,   &
#endif
                                    cvmix_io_close

  Implicit None

!EOP
!BOC

  real(cvmix_r8), parameter :: third = cvmix_one / real(3,cvmix_r8)

  ! CVMix datatypes
  type(cvmix_data_type)       :: CVmix_vars1, CVmix_vars4, CVmix_vars5

  real(cvmix_r8), dimension(:),   allocatable, target :: Mdiff, Tdiff, Sdiff
  real(cvmix_r8), dimension(:),   allocatable, target :: zt, zw_iface,        &
                                                         Ri_bulk, Ri_bulk2
  real(cvmix_r8), dimension(:),   allocatable, target :: w_m, w_s, zeta
  real(cvmix_r8), dimension(:,:), allocatable, target :: TwoDArray
  real(cvmix_r8), dimension(:),   allocatable, target :: buoyancy, shear_sqr, &
                                                         delta_vel_sqr,       &
                                                         buoy_freq_iface
  real(cvmix_r8), dimension(:,:), allocatable, target :: hor_vel
  real(cvmix_r8), dimension(2)                        :: ref_vel
  real(cvmix_r8), dimension(4) :: shape_coeffs
  integer :: i, fid, kt, kw, nlev1, nlev3, nlev4, max_nlev4, OBL_levid4, nlev5
  real(cvmix_r8) :: hmix1, hmix5, ri_crit, layer_thick1, layer_thick4,        &
                    layer_thick5, OBL_depth4, OBL_depth5, N, Nsqr
  real(cvmix_r8) :: kOBL_depth, Bslope, Vslope
  real(cvmix_r8) :: sigma6, OBL_depth6, surf_buoy_force6, surf_fric_vel6,     &
                    vonkarman6, tao, rho0, grav, alpha, Qnonpen, Cp0,         &
                    w_m6, w_s6, wm6_true, ws6_true
  character(len=cvmix_strlen) :: interp_type_t1, interp_type_t4, interp_type_t5
  character(len=cvmix_strlen) :: filename
  ! True => run specified test
  logical :: ltest1, ltest2, ltest3, ltest4, ltest5, ltest6
  logical :: lnoDGat1 ! True => G'(1) = 0 (in test 4)

  namelist/kpp_col1_nml/ltest1, nlev1, layer_thick1, interp_type_t1, hmix1,   &
                        ri_crit
  namelist/kpp_col2_nml/ltest2
  namelist/kpp_col3_nml/ltest3, nlev3
  namelist/kpp_col4_nml/ltest4, interp_type_t4, OBL_levid4, lnoDGat1
  namelist/kpp_col5_nml/ltest5, nlev5, layer_thick5, hmix5, interp_type_t5
  namelist/kpp_col6_nml/ltest6, vonkarman6, tao, rho0, grav, alpha, Qnonpen,  &
                        Cp0, OBL_depth6

  ! Read namelists

  ! Defaults for test 1
  ltest1         = .false.
  nlev1          = 4
  layer_thick1   =  real(10,  cvmix_r8)
  hmix1          = -real(15,  cvmix_r8)
  ri_crit        =  0.3_cvmix_r8
  interp_type_t1 = 'quadratic'

  ! Defaults for test 2
  ltest2 = .false.

  ! Defaults for test 3
  ltest3 = .false.
  nlev3  = 220

  ! Defaults for test 4
  ltest4         = .false.
  OBL_levid4     = 3
  interp_type_t4 = 'quadratic'
  lnoDGat1       = .true.

  ! Defaults for test 5
  ltest5 = .false.
  nlev5 = 10
  layer_thick5   = real(5,  cvmix_r8)
  hmix5          = real(17, cvmix_r8)
  interp_type_t5 = "linear"

  ! Defaults for test 6
  ltest6     = .false.
  vonkarman6 = 0.4_cvmix_r8
  tao        = 0.2_cvmix_r8
  grav       = 9.8_cvmix_r8
  alpha      = 2.5e-4_cvmix_r8
  rho0       =  real(1035, cvmix_r8)
  Qnonpen    = -real(100,  cvmix_r8)
  Cp0        =  real(3992, cvmix_r8)
  OBL_depth6 =  real(6000, cvmix_r8)

  read(*, nml=kpp_col1_nml)
  read(*, nml=kpp_col2_nml)
  read(*, nml=kpp_col3_nml)
  read(*, nml=kpp_col4_nml)
  read(*, nml=kpp_col5_nml)
  read(*, nml=kpp_col6_nml)

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
    call cvmix_put(CVmix_vars1, 'ocn_depth', layer_thick1*real(nlev1,cvmix_r8))

    ! Set up vertical levels (centers and interfaces) and compute bulk
    ! Richardson number
    allocate(zt(nlev1), zw_iface(nlev1+1), Ri_bulk(nlev1))
    do kw=1,nlev1+1
      zw_iface(kw) = -layer_thick1*real(kw-1, cvmix_r8)
    end do
    do kt=1,nlev1
      zt(kt) = 0.5_cvmix_r8*(zw_iface(kt)+zw_iface(kt+1))
      if (zw_iface(kt+1).gt.hmix1) then
        Ri_bulk(kt) = cvmix_zero
      else
        if (Ri_bulk(kt-1).eq.0) then
          ! Exact integration for average value over first cell with non-zero
          ! Ri_bulk
          Ri_bulk(kt) = 0.25_cvmix_r8 * ri_crit *                             &
                        (zw_iface(kt+1)-hmix1)**2 / layer_thick1
        else
          Ri_bulk(kt) = 0.5_cvmix_r8*ri_crit*(hmix1-zt(kt))
        end if
      end if
    end do

    CVmix_vars1%zt_cntr             => zt(:)
    CVmix_vars1%zw_iface            => zw_iface(:)
    CVmix_vars1%BulkRichardson_cntr => Ri_bulk(:)

    ! Compute OBL depth
    call cvmix_kpp_compute_OBL_depth(CVmix_vars1)

    ! Output to screen and file
    print*, "OBL depth = ", CVmix_vars1%BoundaryLayerDepth
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
    call cvmix_output_write_att(fid, "analytic_OBL_depth", -hmix1 +           &
                                                           real(2,cvmix_r8))
    call cvmix_output_write_att(fid, "computed_OBL_depth",                    &
                                CVmix_vars1%BoundaryLayerDepth)
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

    call cvmix_init_kpp(MatchTechnique='MatchGradient')
    call cvmix_kpp_compute_shape_function_coeffs(cvmix_zero, cvmix_zero,      &
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
    call cvmix_init_kpp(vonkarman=cvmix_one, surf_layer_ext=cvmix_one)
    print*, "Coefficients for computing phi_m and phi_s:"
    print*, "a_m = ", cvmix_get_kpp_real('a_m')
    print*, "c_m = ", cvmix_get_kpp_real('c_m')
    print*, "a_s = ", cvmix_get_kpp_real('a_s')
    print*, "c_s = ", cvmix_get_kpp_real('c_s')
    allocate(w_m(nlev3+1), w_s(nlev3+1), zeta(nlev3+1))
    ! Note: zeta = sigma*OBL_depth/MoninObukhov constant
    !       zeta < 0 => unstable flow
    !       zeta > 0 => stable flow
    zeta(1) = -real(2,cvmix_r8)
    do kw=2, nlev3+1
      zeta(kw) = zeta(kw-1) + 2.2_cvmix_r8/real(nlev3,cvmix_r8)
    end do
    ! Typically the first argument of compute_turbulent_scales is sigma, and then
    ! the routine calculates zeta based on the next three parameters. Setting
    ! OBL_depth = surf_buoy_force = surf_fric_vel = 1 (with von Karman = 1 as
    ! well) => sigma = zeta
    call cvmix_kpp_compute_turbulent_scales(zeta, cvmix_one, cvmix_one,       &
                                            cvmix_one, w_m=w_m, w_s=w_s)

    allocate(TwoDArray(nlev3+1,3))
    TwoDArray(:,1) = zeta
    TwoDArray(:,2) = cvmix_one/w_m ! phi_m
    TwoDArray(:,3) = cvmix_one/w_s ! phi_s
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
    deallocate(TwoDArray)
    deallocate(zeta, w_m, w_s)
  endif ! ltest3

  ! Test 4: Compute boundary layer diffusivity
  !         1) Boundary layer between top interface and cell center
  !         2) Boundary layer between cell center and bottom interface
  !         3,4) Same as above, without enhanced diffusivity
  if (ltest4) then
    print*, ""
    print*, "Test 4: Computing Diffusivity in boundary layer"
    print*, "        (2 cases - boundary layer above and below cell center)"
    print*, "        (Both cases run with and without enhanced diffusivity)"
    print*, "----------"

    nlev4 = 5
    max_nlev4 = 10
    if (OBL_levid4.gt.nlev4) &
      OBL_levid4 = nlev4
    layer_thick4 = real(5,cvmix_r8)

    ! Set up vertical levels (centers and interfaces) and compute bulk
    ! Richardson number
    allocate(zt(max_nlev4), zw_iface(max_nlev4+1))
    zt = cvmix_zero
    zw_iface = cvmix_zero
    do kw=1,nlev4+1
      zw_iface(kw) = -layer_thick4*real(kw-1, cvmix_r8)
    end do
    do kt=1,nlev4
      zt(kt) = 0.5_cvmix_r8*(zw_iface(kt)+zw_iface(kt+1))
    end do
    CVmix_vars4%zt_cntr  => zt(:)
    CVmix_vars4%zw_iface => zw_iface(:)

    ! Set up diffusivities
    allocate(Mdiff(max_nlev4+1), Tdiff(max_nlev4+1), Sdiff(max_nlev4+1))
    CVmix_vars4%Mdiff_iface => Mdiff
    CVmix_vars4%Tdiff_iface => Tdiff
    CVmix_vars4%Sdiff_iface => Sdiff

    ! Set physical properties of column for test 4
    call cvmix_put(CVmix_vars4, 'nlev',      nlev4)
    call cvmix_put(CVmix_vars4, 'max_nlev',  max_nlev4)
    call cvmix_put(CVmix_vars4, 'ocn_depth', layer_thick4*real(nlev4,cvmix_r8))
    call cvmix_put(CVmix_vars4, 'surf_fric', cvmix_one)
    call cvmix_put(CVmix_vars4, 'surf_buoy', real(100, cvmix_r8))
    call cvmix_put(CVmix_vars4, 'Coriolis',  1e-4_cvmix_r8)

    do i=1,2
      ! Test 4a/c: Boundary layer above center of level containing it
      Tdiff    = cvmix_zero
      Tdiff(1) = cvmix_one
      Tdiff(2) = real(10,cvmix_r8)
      Tdiff(3) = real(5,cvmix_r8)
      Mdiff = Tdiff
      Sdiff = Tdiff

      OBL_depth4 = abs(zt(OBL_levid4))-layer_thick4/real(4,cvmix_r8)
      call cvmix_put(CVmix_vars4, 'OBL_depth', OBL_depth4)
      call cvmix_put(CVmix_vars4, 'kOBL_depth',                               &
           cvmix_kpp_compute_kOBL_depth(zw_iface, zt, OBL_depth4))

      print*, "OBL_depth = ", CVmix_vars4%BoundaryLayerDepth
      print*, "kOBL_depth = ", CVmix_vars4%kOBL_depth

      call cvmix_init_kpp(ri_crit=ri_crit, vonkarman=0.4_cvmix_r8,            &
                          interp_type2=interp_type_t4, lnoDGat1=lnoDGat1,     &
                          lenhanced_diff=(i.eq.1))
      call cvmix_coeffs_kpp(CVmix_vars4)

      print*, "Height and Diffusivity throughout column: "
      do kw=1,nlev4+1
        write(*,"(1X, F6.2, 1X, F8.5)") zw_iface(kw), Mdiff(kw)
      end do
      print*, ""

      if (i.eq.1) then
        filename="test4a"
      else
        filename="test4c"
      endif
#ifdef _NETCDF
      call cvmix_io_open(fid, trim(filename) // ".nc", "nc")
#else
      call cvmix_io_open(fid, trim(filename) // ".out", "ascii")
#endif

      call cvmix_output_write(fid, CVmix_vars4, (/"zt   ", "zw   ", "Mdiff",  &
                                                  "Tdiff", "Sdiff"/))
#ifdef _NETCDF
      call cvmix_output_write_att(fid, "interp_type2", interp_type_t4)
      call cvmix_output_write_att(fid, "OBL_depth",                           &
                                  CVmix_vars4%BoundaryLayerDepth)
#endif

      call cvmix_io_close(fid)

    ! Test 4b/d: Boundary layer below center of level containing it
      Tdiff    = cvmix_zero
      Tdiff(1) = cvmix_one
      Tdiff(2) = real(10,cvmix_r8)
      Tdiff(3) = real(5,cvmix_r8)
      Mdiff = Tdiff
      Sdiff = Tdiff

      OBL_depth4 = abs(zt(OBL_levid4))+layer_thick4/real(4,cvmix_r8)
      call cvmix_put(CVmix_vars4, 'OBL_depth', OBL_depth4)
      call cvmix_put(CVmix_vars4, 'kOBL_depth',                               &
           cvmix_kpp_compute_kOBL_depth(zw_iface, zt, OBL_depth4))

      print*, "OBL_depth = ", CVmix_vars4%BoundaryLayerDepth
      print*, "kOBL_depth = ", CVmix_vars4%kOBL_depth

      call cvmix_init_kpp(ri_crit=ri_crit, vonkarman=0.4_cvmix_r8,            &
                          interp_type2=interp_type_t4, lnoDGat1=lnoDGat1,     &
                          lenhanced_diff=(i.eq.1))
      call cvmix_coeffs_kpp(CVmix_vars4)

      print*, "Height and Diffusivity throughout column: "
      do kw=1,nlev4+1
        write(*,"(1X, F6.2, 1X, F8.5)") zw_iface(kw), Mdiff(kw)
      end do

      if (i.eq.1) then
        print*, ""
        filename="test4b"
      else
        filename="test4d"
      endif
#ifdef _NETCDF
      call cvmix_io_open(fid, trim(filename) // ".nc", "nc")
#else
      call cvmix_io_open(fid, trim(filename) // ".out", "ascii")
#endif

      call cvmix_output_write(fid, CVmix_vars4, (/"zt   ", "zw   ", "Mdiff",  &
                                                  "Tdiff", "Sdiff"/))
#ifdef _NETCDF
      call cvmix_output_write_att(fid, "interp_type2", interp_type_t4)
      call cvmix_output_write_att(fid, "OBL_depth",                           &
                                  CVmix_vars4%BoundaryLayerDepth)
#endif

      call cvmix_io_close(fid)

    end do

    deallocate(zt, zw_iface)
    deallocate(Mdiff, Tdiff, Sdiff)
  end if ! ltest4

  ! Test 5: Recreate figure C1 from LMD94
  if (ltest5) then
    print*, ""
    print*, "Test 5: Computing Bulk Richardson number"
    print*, "----------"

    ! using linear interpolation, averaging Nsqr, and setting Cv = 1.5  to
    ! match LMD result
    call cvmix_init_kpp(Cv = 1.5_cvmix_r8, interp_type=interp_type_t5)

    ! Set up vertical levels (centers and interfaces) and compute bulk
    ! Richardson number
    allocate(zt(nlev5), zw_iface(nlev5+1))
    do kw=1,nlev5+1
      zw_iface(kw) = -layer_thick5*real(kw-1, cvmix_r8)
    end do
    do kt=1,nlev5
      zt(kt) = 0.5_cvmix_r8 * (zw_iface(kt)+zw_iface(kt+1))
    end do

    ! Compute Br-B(d), |Vr-V(d)|^2, and Vt^2
    allocate(buoyancy(nlev5), delta_vel_sqr(nlev5), hor_vel(nlev5,2),         &
             shear_sqr(nlev5), w_s(nlev5), Ri_bulk(nlev5), Ri_bulk2(nlev5),   &
             buoy_freq_iface(nlev5+1))

    ref_vel(1) = 0.1_cvmix_r8
    ref_vel(2) = cvmix_zero
    N          = 0.01_cvmix_r8
    Nsqr       = N*N
    Bslope     = -Nsqr
    Vslope     = -0.1_cvmix_r8 / (real(nlev5,cvmix_r8)*layer_thick5-hmix5)
    do kt=1,nlev5
      if ((zt(kt).ge.-hmix5).or.(kt.eq.1)) then
        buoyancy(kt)  = Nsqr
        hor_vel(kt,1) = 0.1_cvmix_r8
        buoy_freq_iface(kt) = cvmix_zero
      else
        if (zw_iface(kt).ge.-hmix5) then
          ! derivatives of buoyancy and horizontal velocity component are
          ! discontinuous in this layer (no change -> non-zero linear change)
          ! so we compute area-average of analytic function over layer
          buoyancy(kt)  = Bslope*(-zw_iface(kt+1)-real(hmix5,cvmix_r8))**2 /  &
                          (real(2,cvmix_r8)*layer_thick5) + Nsqr
          hor_vel(kt,1) = Vslope*(-zw_iface(kt+1)-real(hmix5,cvmix_r8))**2 /  &
                          (real(2,cvmix_r8)*layer_thick5) + 0.1_cvmix_r8
        else
          buoyancy(kt)  = Nsqr+Bslope*(-zt(kt)-real(hmix5,cvmix_r8))
          hor_vel(kt,1) = 0.1_cvmix_r8 + Vslope *                             &
                          (-zt(kt)-real(hmix5,cvmix_r8))
        end if
        buoy_freq_iface(kt) = sqrt(-(buoyancy(kt)-buoyancy(kt-1)) /           &
                                    layer_thick5)
      end if
      ! Compute w_s with zeta=0 per LMD page 393
      ! => w_s = von Karman * surf_fric_vel = 0.4*0.01 = 4e-3
      call cvmix_kpp_compute_turbulent_scales(cvmix_zero, -zt(kt),            &
                                              buoyancy(1), 0.01_cvmix_r8,     &
                                              w_s=w_s(kt))
      hor_vel(kt,2) = cvmix_zero
      delta_vel_sqr(kt) = (ref_vel(1)-hor_vel(kt,1))**2 +                     &
                          (ref_vel(2)-hor_vel(kt,2))**2
    end do
    buoy_freq_iface(nlev5+1) = N
!   MNL: test both uses of compute_bulk_Richardson
    Ri_bulk = cvmix_kpp_compute_bulk_Richardson(zt, (buoyancy(1)-buoyancy),   &
                                             delta_vel_sqr,                   &
                                             Nsqr_iface = buoy_freq_iface**2, &
                                             ws_cntr = w_s)

    shear_sqr = cvmix_kpp_compute_unresolved_shear(zt, w_s,                   &
                                             Nsqr_iface = buoy_freq_iface**2)
    ! Note that Vt_shear_sqr is the fourth argument in compute_bulk_Richardson
    ! so it does not need to declared explicitly (even though it is optional)
    Ri_bulk2  = cvmix_kpp_compute_bulk_Richardson(zt, (buoyancy(1)-buoyancy), &
                                                  delta_vel_sqr, shear_sqr)
    call cvmix_kpp_compute_OBL_depth(Ri_bulk, zw_iface, OBL_depth5,           &
                                     kOBL_depth, zt)
    do kt=1,nlev5
      if (abs(Ri_bulk(kt)-Ri_bulk2(kt)).gt.1e-12_cvmix_r8) then
        print*, "WARNING: two Ri_bulk computations did not match!"
        print*, zt(kt), Ri_bulk(kt), Ri_bulk2(kt)
      else
        print*, zt(kt), Ri_bulk(kt)
      end if
    end do
    print*, "OBL has depth of ", OBL_depth5
#ifdef _NETCDF
    print*, "Done! Data is stored in test5.nc, run plot_bulk_Rich.ncl ",      &
            "to see output."
#else
    print*, "Done! Data is stored in test5.out, run plot_bulk_Rich.ncl ",     &
            "to see output."
#endif

    CVmix_vars5%nlev                =  nlev5
    CVmix_vars5%BoundaryLayerDepth  =  OBL_depth5
    CVmix_vars5%zt_cntr             => zt
    CVmix_vars5%BulkRichardson_cntr => Ri_bulk
    CVmix_vars5%Vx_cntr             => hor_vel(:,1)
#ifdef _NETCDF
    call cvmix_io_open(fid, "test5.nc", "nc")
#else
    call cvmix_io_open(fid, "test5.out", "ascii")
#endif
    call cvmix_output_write(fid, CVmix_vars5, (/"zt      ",                   &
                                                "Ri_bulk ",                   &
                                                "Vx      ",                   &
                                                "buoyancy"/), buoyancy)
#ifdef _NETCDF
    call cvmix_output_write_att(fid, "OBL_depth",                             &
                                CVmix_vars5%BoundaryLayerDepth)
    call cvmix_output_write_att(fid, "longname", "ocean height (cell center)",&
                                "zt")
    call cvmix_output_write_att(fid, "units", "m", "zt")
    call cvmix_output_write_att(fid, "longname", "horizontal velocity", "U")
    call cvmix_output_write_att(fid, "units", "m/s", "U")
    call cvmix_output_write_att(fid, "units", "m/s^2", "buoyancy")
    call cvmix_output_write_att(fid, "longname", "Bulk Richardson number",    &
                                "BulkRichardson")
    call cvmix_output_write_att(fid, "units", "unitless", "BulkRichardson")
#endif
    call cvmix_io_close(fid)

    deallocate(zt, zw_iface)
    deallocate(buoyancy, delta_vel_sqr, hor_vel, shear_sqr, w_s, Ri_bulk,     &
               Ri_bulk2, buoy_freq_iface)
  end if ! ltest5

  ! Test 6: Recreate figure C1 from LMD94
  if (ltest6) then
    print*, ""
    print*, "Test 6: 2 simple tests for velocity scale"
    print*, "----------"

    call cvmix_init_kpp(vonkarman=vonkarman6)
    sigma6 = 0.1_cvmix_r8

    surf_buoy_force6 = cvmix_zero
    surf_fric_vel6   = sqrt(tao/rho0)
    print*, "6a: Bf = 0 m^2/s^3 and u* = sqrt(tao/rho0)"
    print*, "                          =", surf_fric_vel6
    wm6_true = cvmix_get_kpp_real("vonkarman")*surf_fric_vel6
    call cvmix_kpp_compute_turbulent_scales(sigma6, OBL_depth6,               &
                                            surf_buoy_force6, surf_fric_vel6, &
                                            w_m = w_m6, w_s = w_s6)

    print*, "    => w_m = w_s ~= vonkarman*u*"
    print*, "                  =", wm6_true
    print*, "w_m = ", w_m6
    print*, "w_s = ", w_s6
    print*, ""

    surf_buoy_force6 = grav*alpha*Qnonpen / (rho0*Cp0)
    surf_fric_vel6   = cvmix_zero
    print*, "6b: u* = 0 m/s and Bf = (grav*alpha/(rho0*Cp0))*Qnonpen"
    print*, "                      =", surf_buoy_force6
    wm6_true = cvmix_get_kpp_real("vonkarman")*                               &
               ((-surf_buoy_force6*OBL_depth6)**third)*                       &
   ((cvmix_get_kpp_real("c_m")*sigma6*cvmix_get_kpp_real("vonkarman"))**third)
    ws6_true = cvmix_get_kpp_real("vonkarman")*                               &
               ((-surf_buoy_force6*OBL_depth6)**third)*                       &
   ((cvmix_get_kpp_real("c_s")*sigma6*cvmix_get_kpp_real("vonkarman"))**third)
    call cvmix_kpp_compute_turbulent_scales(sigma6, OBL_depth6,               &
                                            surf_buoy_force6, surf_fric_vel6, &
                                            w_m = w_m6, w_s = w_s6)

    print*, "    => w_m = vonkarman * (-Bf * OBL_depth)^1/3 *",               &
                          " (c_m*0.1*vonkarman)^1/3 m/s"
    print*, "           = ", wm6_true
    print*, "    => w_s = vonkarman * (-Bf * OBL_depth)^1/3 *",               &
                          " (c_s*0.1*vonkarman)^1/3 m/s"
    print*, "           = ", ws6_true
    print*, "w_m = ", w_m6
    print*, "w_s = ", w_s6
    print*, ""

  end if ! ltest6

!EOC

End Subroutine cvmix_kpp_driver
