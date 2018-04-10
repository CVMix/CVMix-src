!BOP
!\newpage
! !ROUTINE: cvmix_kpp_langmuir_driver
! Qing Li, 20160128

! !DESCRIPTION: A routine to test the KPP module with Langmuir mixing.
!\\
!\\

! !INTERFACE:

subroutine cvmix_kpp_langmuir_driver()

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
   use cvmix_io,              only : cvmix_io_open,                            &
                                     cvmix_input_read,                         &
                                     cvmix_output_write,                       &
                                     cvmix_io_close

   implicit none

!EOP
!BOC

   real(cvmix_r8), parameter :: epssfc = 0.1_cvmix_r8,                         &
                                eps = 1.0e-10_cvmix_r8

   ! CVMix datatypes
   type(cvmix_data_type)       :: CVmix_vars

   real(cvmix_r8), dimension(:),   allocatable, target :: zt, zw_iface,        &
                                                          dz, Ri_bulk
   real(cvmix_r8), dimension(:),   allocatable, target :: w_s
   real(cvmix_r8), dimension(:),   allocatable, target :: bbb, uuu, vvv,       &
                                                          delta_vel_sqr,       &
                                                          buoy_freq_iface,     &
                                                          delta_buoy
   real(cvmix_r8), dimension(:),   allocatable :: ustar_array, u10_array,      &
                                                  b0_array, fcor_array
   real(cvmix_r8), dimension(:,:), allocatable :: bbb_array, uuu_array,        &
                                                  vvv_array
   real(cvmix_r8), dimension(:,:), allocatable :: hbl_kpp

   integer :: ncid_in, ic
   integer :: rfid, kt, kl, kref, ktmp
   integer :: nlev, max_nlev, ncase, jerlov_water_type
   real(cvmix_r8) :: uref, vref, bref, ri_crit, &
                     OBL_depth, surfthick
   real(cvmix_r8) :: kOBL_depth, sigma, bfsfc, ustar, lamult, stable, &
                     fcor, b0, b0sol, absorb_frac, lon, lat, dayofyear, &
                     u10, efactor, init_hbl
   character(len=cvmix_strlen) :: interp_type, MatchTechnique
   character(len=cvmix_strlen) :: infile, outfile, infileEF
   character(len=cvmix_strlen), dimension(:), allocatable :: langmuir_opt
   logical :: llangmuir_efactor, llangmuir_entr, lnoDGat1, luse_efactor_model, l_debug
   type(cvmix_global_params_type) :: CVmix_params

   namelist/langmuir_col_nml/nlev, max_nlev, &
                             interp_type, ri_crit, init_hbl, &
                             llangmuir_efactor, llangmuir_entr, lamult, &
                             jerlov_water_type, &
                             b0, b0sol, ustar, lon, lat, infile, &
                             lnoDGat1, outfile, dayofyear, infileEF

   ! Read namelists

   ! Defaults for Langmuir test
   nlev          = 160
   max_nlev      = 160
   ncase         = 63
   ri_crit       =  0.3_cvmix_r8
   lamult        = cvmix_one
   lon           = 200.0_cvmix_r8
   lat           = 30.0_cvmix_r8
   dayofyear     = 35_cvmix_r8
   init_hbl      = 42_cvmix_r8
   jerlov_water_type = 3
   b0            = -1.190476e-06_cvmix_r8
   b0sol         = cvmix_zero
   interp_type   = 'quadratic'
   infile        = '../../inputdata/langmuir/les_profile_data_20180130.nc'
   outfile       = 'langmuir_test.out'
   lnoDGat1      = .true.
   infileEF      = 'efactor_ww3a_160410.dat'
   MatchTechnique = 'MatchGradient'
   l_debug       = .true.

   ! read(*, nml=langmuir_col_nml)

   ! Allocate variables
   allocate(zt(nlev), zw_iface(nlev+1), dz(nlev))
   allocate(bbb(nlev), uuu(nlev), vvv(nlev), delta_vel_sqr(nlev), w_s(nlev),  &
            Ri_bulk(nlev), buoy_freq_iface(nlev+1), delta_buoy(nlev))
   allocate(ustar_array(ncase), u10_array(ncase), b0_array(ncase),            &
            fcor_array(ncase))
   allocate(bbb_array(nlev,ncase), uuu_array(nlev,ncase),                     &
            vvv_array(nlev,ncase))
   allocate(hbl_kpp(ncase,3))
   allocate(langmuir_opt(3))
   langmuir_opt = (/'off ', 'lw16', 'lf17'/)

   print*, "Test starts"
   print*, "-----------"

!-----------------------------------------------------------------------
! Read input profile data from LES
!-----------------------------------------------------------------------

   ! open netCDF file
   call cvmix_io_open(ncid_in, trim(infile), 'nc', read_only=.true.)
   ! read depth
   call cvmix_input_read(ncid_in, 'z_u', zt)
   call cvmix_input_read(ncid_in, 'z_w', zw_iface(2:nlev+1))
   zw_iface(1) = cvmix_zero
   do kl = 1,nlev
      dz(kl) = zw_iface(kl)-zw_iface(kl+1)
   end do
   ! read surface forcing
   call cvmix_input_read(ncid_in, 'utau',  ustar_array)
   call cvmix_input_read(ncid_in, 'u10',   u10_array)
   call cvmix_input_read(ncid_in, 'bfsfc', b0_array)
   call cvmix_input_read(ncid_in, 'fcor',  fcor_array)
   ! read buoyancy and horizontal velocity profiles
   call cvmix_input_read(ncid_in, 'bxym', bbb_array)
   call cvmix_input_read(ncid_in, 'uxym', uuu_array)
   call cvmix_input_read(ncid_in, 'vxym', vvv_array)
   ! close netCDF file
   call cvmix_io_close(rfid)

!-----------------------------------------------------------------------
! Initialize KPP
!-----------------------------------------------------------------------

   call cvmix_init_kpp(ri_crit=ri_crit, interp_type=interp_type,              &
                      llangmuirEF=llangmuir_efactor,                          &
                      lenhanced_entr=llangmuir_entr,                          &
                      MatchTechnique=MatchTechnique,                          &
                      lnoDGat1=lnoDGat1)
   call cvmix_put(CVmix_vars, 'nlev', nlev)
   call cvmix_put(CVmix_vars, 'max_nlev', max_nlev)
   call cvmix_put(CVmix_vars, 'ocn_depth', -zw_iface(nlev+1))
   call cvmix_put(CVmix_vars, 'Coriolis', fcor)

   CVmix_vars%zt_cntr   => zt(:)
   CVmix_vars%zw_iface  => zw_iface(:)

!-----------------------------------------------------------------------
! Loop over cases
!-----------------------------------------------------------------------

   do ic = 1,ncase
      b0    = b0_array(ic)
      ustar = ustar_array(ic)
      fcor  = fcor_array(ic)
      bbb   = bbb_array(:,ic)
      uuu   = uuu_array(:,ic)
      vvv   = vvv_array(:,ic)

      ! surface buoyancy flux
      call sw_absorb_frac(zt(1), jerlov_water_type, absorb_frac)
      bfsfc = b0 + b0sol*(cvmix_one - absorb_frac)
      stable = merge(cvmix_one, cvmix_zero, bfsfc >= cvmix_zero)
      bfsfc = bfsfc + stable * eps ! ensures bfsfc never =0

      ! get the enhancement factor
      luse_efactor_model = .false.
      select case (trim(langmuir_opt(1)))
         case ('read')
            efactor = cvmix_kpp_efactor_read(infileEF, &
                      lon, lat, dayofyear)
         case ('model')
            efactor = cvmix_kpp_efactor_model(u10, ustar, &
                      init_hbl, CVmix_params)
            luse_efactor_model = .true.
         case ('const')
            efactor = lamult
         case DEFAULT
            efactor = cvmix_one
      end select
      if (l_debug) then
         print*,'efactor0 = ', efactor
      end if

      ! compute bulk Richardson number
      uref = uuu(1)
      vref = vvv(1)
      bref = bbb(1)
      delta_vel_sqr(1) = cvmix_zero
      delta_buoy(1) = cvmix_zero
      buoy_freq_iface(1) = cvmix_zero
      do kl=2,nlev
         ! determin which layer contains surface layer (zt < 0)
         surfthick = epssfc*zt(kl)
         kref = kl
         do ktmp = 1,kl
           if (zw_iface(ktmp+1) .le. surfthick) then
             kref = ktmp
             exit
           end if
         end do

         ! compute uref and vref (layer thickness-weighted average of
         ! layers 1...kref)
         if (kref > 1) then
           uref = uuu(kref)*(-surfthick+zw_iface(kref))
           vref = vvv(kref)*(-surfthick+zw_iface(kref))
           bref = bbb(kref)*(-surfthick+zw_iface(kref))
           do ktmp = 1,kref-1
             uref = uref + dz(ktmp)*uuu(ktmp)
             vref = vref + dz(ktmp)*vvv(ktmp)
             bref = bref + dz(ktmp)*bbb(ktmp)
           end do
           uref = -uref / surfthick
           vref = -vref / surfthick
           bref = -bref / surfthick
         else
           uref = uuu(1)
           vref = vvv(1)
           bref = bbb(1)
         end if

         ! calculate the absorbed bouyancy flux
         call sw_absorb_frac(zt(kl),jerlov_water_type,absorb_frac)
         bfsfc = b0 + b0sol*(cvmix_one - absorb_frac)
         stable = merge(cvmix_one, cvmix_zero, bfsfc >= cvmix_zero)
         bfsfc = bfsfc + stable * eps ! ensures bfsfc never =0

         ! buoyancy differences
         delta_buoy(kl) = bref-bbb(kl)

         ! velocity shear
         delta_vel_sqr(kl) = (uref - uuu(kl))**2 + (vref - vvv(kl))**2

         ! buoyancy frequency
         buoy_freq_iface(kl) = sqrt(2.0_cvmix_r8 *                     &
                               (bbb(kl-1)-bbb(kl))/(dz(kl)+dz(kl-1)))

         ! compute velocity scales at sigma, for hbl = -zt(kl)
         sigma = epssfc
         call cvmix_kpp_compute_turbulent_scales(sigma, -zt(kl),       &
                 bfsfc, ustar, langmuir_Efactor=efactor, w_s=w_s(kl))
         !print*, "kl = ", kl, "ws = ",w_s(kl)
      end do

      ! compute bulk Richardson number
      Ri_bulk = cvmix_kpp_compute_bulk_Richardson(zt, delta_buoy,      &
                                                  delta_vel_sqr,       &
                                     Nsqr_iface = buoy_freq_iface**2,  &
                                        ws_cntr = w_s)
      ! compute boundary layer depth
      call cvmix_kpp_compute_OBL_depth(Ri_bulk, zw_iface, OBL_depth,   &
                                     kOBL_depth, zt)

      print*, "Bulk Richardson number "
      do kt=1,nlev
        print*, zt(kt), Ri_bulk(kt)
      end do
      print*, " "
      print*, "OBL has depth of ", OBL_depth
      print*, " "

      print*, "Done! Data is stored in ", trim(outfile)
   end do ! loop over cases
   ! call cvmix_io_open(fid, trim(outfile), "ascii")
   ! do kt=1,nlev+1
   !   write(fid,100) zw_iface(kt), uw(kt),vw(kt),wb(kt)
   ! end do
   ! call cvmix_io_close(fid)

   deallocate(zt, zw_iface)
   deallocate(uuu, vvv, bbb, delta_vel_sqr, w_s, Ri_bulk, buoy_freq_iface)
   ! deallocate(Mdiff, Tdiff, Sdiff, Ghat)
   ! deallocate(uw,vw,wb)
!EOC

end subroutine cvmix_kpp_langmuir_driver

!-----------------------------------------------------------------------

subroutine sw_absorb_frac(depth, jerlov_water_type, sw_absorb_fraction)

! !DESCRIPTION:
!  Computes fraction of solar short-wave flux penetrating to
!  specified depth due to exponential decay in Jerlov water type.
!  Reference : two band solar absorption model of Simpson and
!     Paulson (1977)
!  Note: below 200m the solar penetration gets set to zero,
!     otherwise the limit for the exponent ($+/- 5678$) needs to be
!     taken care of.
!
! !REVISION HISTORY:
!  Borrowed from POP2


   use cvmix_kinds_and_types, only : cvmix_r8,                         &
                                     cvmix_zero,                       &
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
      depth_cutoff = -200.0

!-----------------------------------------------------------------------
!
!   define Jerlov water properties with rfac, depth1, depth2
!     Jerlov water type :  I       IA      IB      II      III
!     jerlov_water_type :  1       2       3       4       5
!
!-----------------------------------------------------------------------

   real (cvmix_r8), dimension(num_water_types) ::                      &
      rfac   = (/ 0.58, 0.62, 0.67, 0.77, 0.78 /),                     &
      depth1 = (/ 0.35, 0.60, 1.00, 1.50, 1.40 /),                     &
      depth2 = (/ 23.0, 20.0, 17.0, 14.0, 7.90 /)

!-----------------------------------------------------------------------
!
!  compute absorption fraction
!
!-----------------------------------------------------------------------

   if (depth < depth_cutoff) then
      sw_absorb_fraction = cvmix_zero
   else
      sw_absorb_fraction = rfac(jerlov_water_type) *                   &
                 exp(depth/depth1(jerlov_water_type)) +                &
                     (cvmix_one - rfac(jerlov_water_type))*            &
                 exp(depth/depth2(jerlov_water_type))
   endif

end subroutine sw_absorb_frac
