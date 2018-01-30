!BOP
!\newpage
! !ROUTINE: cvmix_kpp_langmuir_driver
! Qing Li, 20160128

! !DESCRIPTION: A routine to test the KPP module with Langmuir mixing.
!\\
!\\

! !INTERFACE:

Subroutine cvmix_kpp_langmuir_driver()

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                 &
                                    cvmix_zero,               &
                                    cvmix_one,                &
                                    cvmix_strlen,             &
                                    cvmix_data_type,          &
                                    cvmix_global_params_type
  use cvmix_kpp,             only : cvmix_init_kpp,                           &
                                    cvmix_get_kpp_real,                       &
                                    cvmix_kpp_compute_OBL_depth,              &
                                    cvmix_kpp_compute_kOBL_depth,             &
                                    cvmix_kpp_compute_bulk_Richardson,        &
                                    cvmix_kpp_compute_unresolved_shear,       &
                                    cvmix_kpp_compute_turbulent_scales,       &
                                    cvmix_kpp_compute_shape_function_coeffs,  &
                                    cvmix_kpp_efactor_read,                   &
                                    cvmix_kpp_efactor_model,                  &
                                    cvmix_coeffs_kpp
  use cvmix_put_get,         only : cvmix_put
  use cvmix_io,              only : cvmix_io_open,            &
                                    cvmix_output_write,       &
                                    cvmix_io_close

  Implicit None

!EOP
!BOC

  real(cvmix_r8), parameter :: epssfc = 0.1_cvmix_r8, &
                               eps = 1.0e-10_cvmix_r8, &
                               epsz = 1.0e-3_cvmix_r8

  ! CVMix datatypes
  type(cvmix_data_type)       :: CVmix_vars

  real(cvmix_r8), dimension(:),   allocatable, target :: Mdiff, Tdiff, &
                                                         Sdiff, Ghat
  real(cvmix_r8), dimension(:),   allocatable, target :: zt, zw_iface, &
                                                         zz, Ri_bulk
  real(cvmix_r8), dimension(:),   allocatable, target :: w_s
  real(cvmix_r8), dimension(:),   allocatable, target :: bbb, uuu, vvv, &
                                                         delta_vel_sqr, &
                                                         buoy_freq_iface, &
                                                         delta_buoy 
  real(cvmix_r8), dimension(:),   allocatable :: uw, vw, wb
  integer :: fid, rfid, kt, kw, kl, kref, ktmp, nlev, max_nlev, &
             jerlov_water_type
  real(cvmix_r8) :: uref, vref, bref, ri_crit, layer_thick, &
                    OBL_depth, surfthick
  real(cvmix_r8) :: kOBL_depth, sigma, bfsfc, ustar, lamult, stable, tmp, &
                    Coriolis, b0, b0sol, absorb_frac, lon, lat, dayofyear, &
                    u10, efactor, init_hbl
  character(len=cvmix_strlen) :: interp_type, infile, outfile, infileEF
  character(len=cvmix_strlen) :: optionEF
  logical :: llangmuir_efactor, lnoDGat1, luse_efactor_model
  type(cvmix_global_params_type) :: CVmix_params

  namelist/langmuir_col_nml/nlev, max_nlev, layer_thick, &
                            interp_type, ri_crit, optionEF, u10, init_hbl, &
                            llangmuir_efactor, lamult, jerlov_water_type, & 
                            b0, b0sol,Coriolis, ustar, lon, lat, infile, &
                            lnoDGat1, outfile, dayofyear, infileEF

  ! Read namelists

  ! Defaults for Langmuir test
  nlev          = 5
  max_nlev      = 10
  layer_thick   =  real(10,  cvmix_r8)
  ri_crit       =  0.3_cvmix_r8
  llangmuir_efactor = .false.
  optionEF      = 'default'
  lamult        = cvmix_one
  lon           = 200.0_cvmix_r8
  lat           = 30.0_cvmix_r8
  dayofyear     = 35_cvmix_r8
  ustar         = 6.1e-3_cvmix_r8
  u10           = 5_cvmix_r8
  init_hbl      = 33_cvmix_r8
  jerlov_water_type = 3
  b0            = -1.190476e-06_cvmix_r8
  b0sol         = 1.190476e-06_cvmix_r8
  Coriolis      = 1e-4_cvmix_r8
  interp_type   = 'quadratic'
  infile        = 'uvbPfl.txt'
  outfile       = 'langmuir_test.out'
  lnoDGat1      = .true.
  infileEF      = 'efactor_ww3a_160410.dat'

  read(*, nml=langmuir_col_nml)

  ! Allocate variables
  allocate(zt(nlev), zw_iface(nlev+1), zz(nlev))
  allocate(bbb(nlev), uuu(nlev), vvv(nlev), delta_vel_sqr(nlev), &
           w_s(nlev), Ri_bulk(nlev), buoy_freq_iface(nlev+1), &
           delta_buoy(nlev))
  allocate(Mdiff(max_nlev+1),Tdiff(max_nlev+1),Sdiff(max_nlev+1), &
           Ghat(max_nlev+1))
  allocate(uw(nlev+1),vw(nlev+1),wb(nlev+1))

  ! Read buoyancy, horizontal velocity
  rfid = 2
  call cvmix_io_open(rfid,trim(infile),'ascii',.true.)
  do kt=1,nlev
    read(rfid,*) zz(kt),uuu(kt),vvv(kt),bbb(kt)
  end do
  call cvmix_io_close(rfid)
  ! Chech vertical levels are consistent with setup
  tmp = zz(1)-zz(2)
  if (abs(tmp-layer_thick) .ge. epsz) then
     write(*,*) tmp, layer_thick
     print*, "Vertical levels not consistent. Stop"
     stop
  end if
  !print*,"b = ", bbb

  print*, "Test starts"
  print*, "----------"

  ! Initialize parameter datatype and set up column
  call cvmix_init_kpp(ri_crit=ri_crit, interp_type=interp_type,&
                      llangmuirEF=llangmuir_efactor, &
                      MatchTechnique='MatchGradient', &
                      lnoDGat1=lnoDGat1)
  call cvmix_put(CVmix_vars, 'nlev', nlev)
  call cvmix_put(CVmix_vars, 'max_nlev', max_nlev)
  call cvmix_put(CVmix_vars, 'ocn_depth', layer_thick*real(nlev,cvmix_r8))
  call cvmix_put(CVmix_vars, 'Coriolis', Coriolis)

  ! Set up vertical levels (centers and interfaces) 
  ! interfaces
  do kw=1,nlev+1
    zw_iface(kw) = -layer_thick*real(kw-1, cvmix_r8)
    !print*, zw_iface(kw)
  end do
  ! centers
  do kt=1,nlev
    zt(kt) = 0.5_cvmix_r8*(zw_iface(kt)+zw_iface(kt+1))
    !print*, zt(kt)
  end do
  CVmix_vars%zt_cntr   => zt(:)
  CVmix_vars%zw_iface  => zw_iface(:)

  ! Surface buoyancy flux
  call sw_absorb_frac(zt(1),jerlov_water_type,absorb_frac)
  bfsfc = b0 + b0sol*(cvmix_one - absorb_frac)
  stable = merge(cvmix_one, cvmix_zero, bfsfc >= cvmix_zero)
  bfsfc = bfsfc + stable * eps ! ensures bfsfc never =0
  !print*,'bfsfc = ', bfsfc
  !print*,'stable = ', stable

  ! Get the enhancement factor
  luse_efactor_model = .false.
  select case (trim(optionEF))
    case ('read')
      efactor = cvmix_kpp_efactor_read(infileEF, lon, lat, dayofyear)
    case ('model')
      efactor = cvmix_kpp_efactor_model(u10, ustar, init_hbl, CVmix_params)
      luse_efactor_model = .true.
    case DEFAULT
      efactor = lamult
  end select
  print*,'efactor0 = ', efactor

  ! Compute bulk Richardson number
  uref = uuu(1)
  vref = vvv(1)
  bref = bbb(1)
  delta_vel_sqr(1) = cvmix_zero
  delta_buoy(1) = cvmix_zero
  buoy_freq_iface(1) = cvmix_zero
  do kl=2,nlev
    ! Determin which layer contains surface layer (zt < 0)
    surfthick = epssfc*zt(kl)
    kref = kl
    do ktmp = 1,kl
      if (zw_iface(ktmp+1) .le. surfthick) then
        kref = ktmp
        exit
      end if
    end do
    !print*, kref

    ! Compute uref and vref (layer thickness-weighted average of
    ! layers 1...kref)
    if (kref > 1) then
      uref = uuu(kref)*(-surfthick+zw_iface(kref))
      vref = vvv(kref)*(-surfthick+zw_iface(kref))
      bref = bbb(kref)*(-surfthick+zw_iface(kref))
      do ktmp = 1,kref-1
        uref = uref + layer_thick*uuu(ktmp)
        vref = vref + layer_thick*vvv(ktmp)
        bref = bref + layer_thick*bbb(ktmp)
      end do
      uref = -uref / surfthick
      vref = -vref / surfthick
      bref = -bref / surfthick
    !  uref = uuu(kref)
    !  vref = vvv(kref)
    !  bref = bbb(kref)
    else
      uref = uuu(1)
      vref = vvv(1)
      bref = bbb(1)
    end if
    !print*,'uref = ',uref,' vref = ',vref, ' bref = ',bref
    !print*,'uref2 = ',uuu(kref),' vref2 = ',vvv(kref), ' bref2 = ',bbb(kref)
    
    ! calculate the absorbed bouyancy flux
    call sw_absorb_frac(zt(kl),jerlov_water_type,absorb_frac)
    bfsfc = b0 + b0sol*(cvmix_one - absorb_frac)
    stable = merge(cvmix_one, cvmix_zero, bfsfc >= cvmix_zero)
    bfsfc = bfsfc + stable * eps ! ensures bfsfc never =0
    !print*,'bfsfc = ', bfsfc

    ! Buoyancy differences
    delta_buoy(kl) = bref-bbb(kl)

    ! Velocity shear
    delta_vel_sqr(kl) = (uref - uuu(kl))**2 + (vref - vvv(kl))**2

    ! Buoyancy frequency
    buoy_freq_iface(kl) = sqrt(-(bbb(kl)-bbb(kl-1)) / layer_thick)

    ! Compute velocity scales at sigma, for hbl = -zt(kl)
    sigma = epssfc
    !sigma = cvmix_zero
    call cvmix_kpp_compute_turbulent_scales(sigma, &
            -zt(kl),                &
            bfsfc,        &
            ustar,          &
            langmuir_Efactor=efactor,     &
            w_s=w_s(kl))
    !print*, "kl = ", kl, "ws = ",w_s(kl)
  end do
  ! Compute bulk Richardson number
  Ri_bulk = cvmix_kpp_compute_bulk_Richardson(zt, delta_buoy,   &
                                              delta_vel_sqr,      &
                                 Nsqr_iface = buoy_freq_iface**2, &
                                    ws_cntr = w_s)
  ! Compute boundary layer depth
  call cvmix_kpp_compute_OBL_depth(Ri_bulk, zw_iface, OBL_depth,  &
                                    kOBL_depth, zt)

  print*, "Bulk Richardson number "
  do kt=1,nlev
    print*, zt(kt), Ri_bulk(kt)
  end do
  print*, " "
  print*, "OBL has depth of ", OBL_depth
  print*, " "

  ! update the absorbed bouyancy flux to OBL_depth
  call sw_absorb_frac(-OBL_depth,jerlov_water_type,absorb_frac)
  bfsfc = b0 + b0sol*(cvmix_one - absorb_frac)
  stable = merge(cvmix_one, cvmix_zero, bfsfc >= cvmix_zero)
  bfsfc = bfsfc + stable * eps ! ensures bfsfc never =0
  !print*,'bfsfc = ', bfsfc

  ! update efactor if luse_efactor_model is true
  if (luse_efactor_model) then
    efactor = cvmix_kpp_efactor_model(u10, ustar, OBL_depth, CVmix_params)
  end if
  print*,"efactor1 = ", efactor

  ! Set up diffusivities
  CVmix_vars%BoundaryLayerDepth  =  OBL_depth
  CVmix_vars%kOBL_depth  =  kOBL_depth
  CVmix_vars%BulkRichardson_cntr => Ri_bulk
  CVmix_vars%SurfaceFriction = ustar
  CVmix_vars%SurfaceBuoyancyForcing = bfsfc
  CVmix_vars%LangmuirEnhancementFactor = efactor
  CVmix_vars%Mdiff_iface => Mdiff
  CVmix_vars%Tdiff_iface => Tdiff
  CVmix_vars%Sdiff_iface => Sdiff
  CVmix_vars%kpp_Tnonlocal_iface => Ghat

  Tdiff = cvmix_zero
  Mdiff = Tdiff
  Sdiff = Tdiff
  Ghat = Tdiff

  ! Compute diffusivities
  call cvmix_coeffs_kpp(CVmix_vars)

  !print*, "Height and Diffusivity throughout column: "
  !print*, "z_face    Mdiff    Tdiff     Ghat"
  !do kt=1,nlev+1
  !  write(*,"(1X, F6.2, 1X, F8.5, 1X, F8.5, 1X, F8.5)") &
  !         zw_iface(kt), Mdiff(kt), Tdiff(kt), Ghat(kt)
  !end do
  !print*, " "
 
  ! Compute fluxes
  uw = cvmix_zero
  vw = cvmix_zero
  wb = cvmix_zero
  uw(1) = ustar**2
  wb(1) = -(b0+b0sol)+stable*eps
  do kt=2,nlev
    uw(kt) = -Mdiff(kt)*(uuu(kt-1)-uuu(kt))/layer_thick
    vw(kt) = -Mdiff(kt)*(vvv(kt-1)-vvv(kt))/layer_thick
    wb(kt) = -Tdiff(kt)*(bbb(kt-1)-bbb(kt))/layer_thick &
             -Ghat(kt)*bfsfc
  end do

  !print*, "Height and vertical turbulent fluxes: "
  !print*, "z_face              uw              vw              wb"
  !do kt=1,nlev+1
  !  write(*,100) zw_iface(kt), uw(kt), vw(kt), wb(kt)
  !end do
  !print*, " "
100 format(1X, F6.2, 1X, E15.7, 1X, E15.7, 1X, E15.7) 

  print*, "Done! Data is stored in ", trim(outfile)

  call cvmix_io_open(fid, trim(outfile), "ascii")
  do kt=1,nlev+1
    write(fid,100) zw_iface(kt), uw(kt),vw(kt),wb(kt)
  end do
  call cvmix_io_close(fid)

  deallocate(zt, zw_iface)
  deallocate(uuu, vvv, bbb, delta_vel_sqr, w_s, Ri_bulk, buoy_freq_iface)
  deallocate(Mdiff, Tdiff, Sdiff, Ghat)
  deallocate(uw,vw,wb)
!EOC

End Subroutine cvmix_kpp_langmuir_driver

subroutine sw_absorb_frac( depth, jerlov_water_type, sw_absorb_fraction )

!  Computes fraction of solar short-wave flux penetrating to
!  specified depth due to exponential decay in Jerlov water type.
!  Reference : two band solar absorption model of Simpson and
!     Paulson (1977)
!  Note: below 200m the solar penetration gets set to zero,
!     otherwise the limit for the exponent ($+/- 5678$) needs to be 
!     taken care of.
!

   use cvmix_kinds_and_types, only : cvmix_r8,                 &
                                    cvmix_zero,               &
                                    cvmix_one
! !INPUT PARAMETERS:

   real (cvmix_r8), intent(in) :: &
      depth     ! vertical depth (m, <0.) for desired sw fraction
      
   integer, intent(in) :: jerlov_water_type ! Jerlov water type

! !OUTPUT PARAMETERS:

   real (cvmix_r8), intent(out) :: &
     sw_absorb_fraction     ! short wave (radiation) fractional decay

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer, parameter :: &
      num_water_types = 5  ! max number of different water types

   real (cvmix_r8), parameter :: &
      depth_cutoff = -200.0_cvmix_r8

!-----------------------------------------------------------------------
!
!   define Jerlov water properties with rfac, depth1, depth2
!     Jerlov water type :  I       IA      IB      II      III
!     jerlov_water_type :  1       2       3       4       5
!
!-----------------------------------------------------------------------

   real (cvmix_r8), dimension(num_water_types) ::                       &
      rfac   = (/ 0.58_cvmix_r8, 0.62_cvmix_r8, 0.67_cvmix_r8, 0.77_cvmix_r8, 0.78_cvmix_r8 /), &
      depth1 = (/ 0.35_cvmix_r8, 0.60_cvmix_r8, 1.00_cvmix_r8, 1.50_cvmix_r8, 1.40_cvmix_r8 /), &
      depth2 = (/ 23.0_cvmix_r8, 20.0_cvmix_r8, 17.0_cvmix_r8, 14.0_cvmix_r8, 7.90_cvmix_r8 /)

!-----------------------------------------------------------------------
!
!  compute absorption fraction
!
!-----------------------------------------------------------------------

   if (depth < depth_cutoff) then
      sw_absorb_fraction = cvmix_zero
   else
      sw_absorb_fraction =     rfac(jerlov_water_type)*            &
                 exp(depth/depth1(jerlov_water_type)) +  &
                     (cvmix_one - rfac(jerlov_water_type))*                &
                 exp(depth/depth2(jerlov_water_type))
   endif

 end subroutine sw_absorb_frac
