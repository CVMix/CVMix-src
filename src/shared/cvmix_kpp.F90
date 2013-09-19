!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module cvmix_kpp

!BOP
!\newpage
! !MODULE: cvmix_kpp
!
! !DESCRIPTION:
!  This module contains routines to initialize the derived types needed for
!  KPP mixing and to set the viscosity and diffusivity coefficients
!  accordingly.
!\\
!\\
!
! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                     &
                                    one,                          &
                                    cvmix_data_type
  use cvmix_put_get, only :         cvmix_put
  use cvmix_math, only :            CVMIX_MATH_INTERP_LINEAR,      &
                                    CVMIX_MATH_INTERP_QUAD,        &
                                    CVMIX_MATH_INTERP_CUBE_SPLINE, &
                                    cvmix_math_poly_interp,        &
                                    cvmix_math_cubic_root_find,    &
                                    cvmix_math_evaluate_cubic

!EOP

  implicit none
  private
  save

!BOP

! !DEFINED_PARAMETERS:
  integer, parameter :: CVMIX_KPP_INTERP_POP = -1

! !PUBLIC MEMBER FUNCTIONS:

  public :: cvmix_init_kpp
  ! Note: cvmix_kpp_compute_OBL_depth would be part of cvmix_coeffs_kpp but
  !       CVMix can not smooth the boundary layer depth or correct the
  !       buoyancy flux term
  public :: cvmix_kpp_compute_OBL_depth
  public :: cvmix_coeffs_kpp
  public :: cvmix_put_kpp
  public :: cvmix_get_kpp_real
  public :: cvmix_kpp_compute_bulk_Richardson
  public :: cvmix_kpp_compute_turbulent_scales
  public :: cvmix_kpp_compute_unresolved_shear
  ! These are public for testing, may end up private later
  public :: cvmix_kpp_compute_shape_function_coeffs
  public :: cvmix_kpp_compute_kOBL_depth

  interface cvmix_coeffs_kpp
    module procedure cvmix_coeffs_kpp_low
    module procedure cvmix_coeffs_kpp_wrap
  end interface cvmix_coeffs_kpp

  interface cvmix_put_kpp
    module procedure cvmix_put_kpp_int
    module procedure cvmix_put_kpp_real
    module procedure cvmix_put_kpp_logical
  end interface cvmix_put_kpp

  interface cvmix_kpp_compute_OBL_depth
    module procedure cvmix_kpp_compute_OBL_depth_low
    module procedure cvmix_kpp_compute_OBL_depth_wrap
  end interface cvmix_kpp_compute_OBL_depth

  interface cvmix_kpp_compute_turbulent_scales
    module procedure cvmix_kpp_compute_turbulent_scales_0d
    module procedure cvmix_kpp_compute_turbulent_scales_1d
  end interface cvmix_kpp_compute_turbulent_scales

! !PUBLIC TYPES:

  ! cvmix_kpp_params_type contains the necessary parameters for KPP mixing
  type, public :: cvmix_kpp_params_type
    private
    real(cvmix_r8) :: Ri_crit        ! Critical Richardson number
                                     ! (OBL_depth = where bulk Ri = Ri_crit)
    real(cvmix_r8) :: vonkarman      ! von Karman constant
    real(cvmix_r8) :: Cstar          ! coefficient for nonlinear transport
    ! For velocity scale function, _m => momentum and _s => scalar (tracer)
    real(cvmix_r8) :: zeta_m         ! parameter for computing vel scale func
    real(cvmix_r8) :: zeta_s         ! parameter for computing vel scale func
    real(cvmix_r8) :: a_m            ! parameter for computing vel scale func
    real(cvmix_r8) :: a_s            ! parameter for computing vel scale func
    real(cvmix_r8) :: c_m            ! parameter for computing vel scale func
    real(cvmix_r8) :: c_s            ! parameter for computing vel scale func
    real(cvmix_r8) :: surf_layer_ext ! nondimensional extent of surface layer
                                     ! (expressed in sigma-coordinates)
    integer        :: interp_type    ! interpolation type used to interpolate
                                     ! bulk Richardson number
    integer        :: interp_type2   ! interpolation type used to interpolate
                                     ! diff and visc at OBL_depth
    logical        :: lEkman         ! True => compute Ekman depth limit
    logical        :: lMonOb         ! True => compute Monin-Obukhov limit
    logical        :: lnoDGat1       ! True => G'(1) = 0 (shape function)
                                     ! False => compute G'(1) as in LMD94
    logical        :: lavg_N_or_Nsqr ! True => N (or Nsqr) at cell center is
                                     !  average of values at interfaces above
                                     !  and below.
                                     ! False => N (or Nsqr) at cell center is
                                     !  set to value at interface below
                                     ! (only used in compute_unresolved_shear)
  end type cvmix_kpp_params_type

!EOP

type(cvmix_kpp_params_type), target :: CVmix_kpp_params_saved

contains

!BOP

! !IROUTINE: cvmix_init_kpp
! !INTERFACE:

  subroutine cvmix_init_kpp(ri_crit, vonkarman, Cstar, zeta_m, zeta_s, a_m,   &
                            a_s, c_m, c_s, surf_layer_ext, interp_type,       &
                            interp_type2, lEkman, lMonOb, lnoDGat1,           &
                            lavg_N_or_Nsqr, CVmix_kpp_params_user)

! !DESCRIPTION:
!  Initialization routine for KPP mixing.
!
! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    real(cvmix_r8),   optional :: ri_crit, &     ! units: unitless
                                  vonkarman, &   ! units: unitless
                                  Cstar, &       ! units: unitless
                                  zeta_m, &      ! units: unitless
                                  zeta_s, &      ! units: unitless
                                  a_m, &         ! units: unitless
                                  a_s, &         ! units: unitless
                                  c_m, &         ! units: unitless
                                  c_s, &         ! units: unitless
                                  surf_layer_ext ! units: unitless
    character(len=*), optional :: interp_type, interp_type2
    logical,          optional :: lEkman, lMonOb, lnoDGat1, lavg_N_or_Nsqr

! !OUTPUT PARAMETERS:
    type(cvmix_kpp_params_type), intent(inout), target, optional ::           &
                                              CVmix_kpp_params_user

!EOP
!BOC

    type(cvmix_kpp_params_type), pointer :: CVmix_kpp_params_out

    CVmix_kpp_params_out => CVmix_kpp_params_saved
    if (present(CVmix_kpp_params_user)) then
      CVmix_kpp_params_out => CVmix_kpp_params_user
    end if

    if (present(ri_crit)) then
      call cvmix_put_kpp('Ri_crit', ri_crit, CVmix_kpp_params_user)
    else
      call cvmix_put_kpp('Ri_crit', 0.3_cvmix_r8, CVmix_kpp_params_user)
    end if

    if (present(vonkarman)) then
      call cvmix_put_kpp('vonkarman', vonkarman, CVmix_kpp_params_user)
    else
      call cvmix_put_kpp('vonkarman', 0.41_cvmix_r8, CVmix_kpp_params_user)
    end if

    if (present(Cstar)) then
      call cvmix_put_kpp('Cstar', Cstar, CVmix_kpp_params_user)
    else
      call cvmix_put_kpp('Cstar', 10.0_cvmix_r8, CVmix_kpp_params_user)
    end if

    if (present(zeta_m)) then
      call cvmix_put_kpp('zeta_m', zeta_m, CVmix_kpp_params_user)
    else
      call cvmix_put_kpp('zeta_m', -0.2_cvmix_r8, CVmix_kpp_params_user)
    end if

    if (present(zeta_s)) then
      call cvmix_put_kpp('zeta_s', zeta_s, CVmix_kpp_params_user)
    else
      call cvmix_put_kpp('zeta_s', -1.0_cvmix_r8, CVmix_kpp_params_user)
    end if

    if (present(a_m)) then
      call cvmix_put_kpp('a_m', a_m, CVmix_kpp_params_user)
    else
      call cvmix_put_kpp('a_m', 1.26_cvmix_r8, CVmix_kpp_params_user)
    end if

    if (present(a_s)) then
      call cvmix_put_kpp('a_s', a_s, CVmix_kpp_params_user)
    else
      call cvmix_put_kpp('a_s', -28.86_cvmix_r8, CVmix_kpp_params_user)
    end if

    if (present(c_m)) then
      call cvmix_put_kpp('c_m', c_m, CVmix_kpp_params_user)
    else
      call cvmix_put_kpp('c_m', 8.38_cvmix_r8, CVmix_kpp_params_user)
    end if

    if (present(c_s)) then
      call cvmix_put_kpp('c_s', c_s, CVmix_kpp_params_user)
    else
      call cvmix_put_kpp('c_s', 98.96_cvmix_r8, CVmix_kpp_params_user)
    end if

    if (present(surf_layer_ext)) then
      call cvmix_put_kpp('surf_layer_ext', surf_layer_ext,                    &
                         CVmix_kpp_params_user)
    else
      call cvmix_put_kpp('surf_layer_ext', 0.1_cvmix_r8, CVmix_kpp_params_user)
    end if

    if (present(interp_type)) then
      select case (trim(interp_type))
        case ('line', 'linear')
          call cvmix_put_kpp(CVmix_kpp_params_out, 'interp_type', &
                             CVMIX_MATH_INTERP_LINEAR)
        case ('quad', 'quadratic')
          call cvmix_put_kpp(CVmix_kpp_params_out, 'interp_type', &
                             CVMIX_MATH_INTERP_QUAD)
        case ('cube', 'cubic', 'cubic_spline', 'cubic spline')
          call cvmix_put_kpp(CVmix_kpp_params_out, 'interp_type', &
                             CVMIX_MATH_INTERP_CUBE_SPLINE)
        case DEFAULT
          print*, "ERROR: ", trim(interp_type), " is not a valid type of ", &
                  "interpolation!"
          stop 1
      end select
    else
      call cvmix_put_kpp(CVmix_kpp_params_out, 'interp_type', &
                         CVMIX_MATH_INTERP_QUAD)
    end if

    if (present(interp_type2)) then
      select case (trim(interp_type2))
        case ('line', 'linear')
          call cvmix_put_kpp(CVmix_kpp_params_out, 'interp_type2', &
                             CVMIX_MATH_INTERP_LINEAR)
        case ('quad', 'quadratic')
          call cvmix_put_kpp(CVmix_kpp_params_out, 'interp_type2', &
                             CVMIX_MATH_INTERP_QUAD)
        case ('cube', 'cubic', 'cubic_spline', 'cubic spline')
          call cvmix_put_kpp(CVmix_kpp_params_out, 'interp_type2', &
                             CVMIX_MATH_INTERP_CUBE_SPLINE)
        case ('POP')
          call cvmix_put_kpp(CVmix_kpp_params_out, 'interp_type2', &
                             CVMIX_KPP_INTERP_POP)
        case DEFAULT
          print*, "ERROR: ", trim(interp_type), " is not a valid type of ", &
                  "interpolation!"
          stop 1
      end select
    else
      call cvmix_put_kpp(CVmix_kpp_params_out, 'interp_type2', &
                         CVMIX_MATH_INTERP_CUBE_SPLINE)
    end if

    if (present(lEkman)) then
      call cvmix_put_kpp(CVmix_kpp_params_out, 'lEkman', lEkman)
    else
      call cvmix_put_kpp(CVmix_kpp_params_out, 'lEkman', .false.)
    end if

    if (present(lMonOb)) then
      call cvmix_put_kpp(CVmix_kpp_params_out, 'lMonOb', lMonOb)
    else
      call cvmix_put_kpp(CVmix_kpp_params_out, 'lMonOb', .false.)
    end if

    if (present(lnoDGat1)) then
      call cvmix_put_kpp(CVmix_kpp_params_out, 'lnoDGat1', lnoDGat1)
    else
      call cvmix_put_kpp(CVmix_kpp_params_out, 'lnoDGat1', .true.)
    end if

    if (present(lavg_N_or_Nsqr)) then
      call cvmix_put_kpp(CVmix_kpp_params_out,'lavg_N_or_Nsqr',lavg_N_or_Nsqr)
    else
      call cvmix_put_kpp(CVmix_kpp_params_out,'lavg_N_or_Nsqr',.false.)
    end if

!EOC

  end subroutine cvmix_init_kpp

!BOP

! !IROUTINE: cvmix_coeffs_kpp_wrap
! !INTERFACE:

  subroutine cvmix_coeffs_kpp_wrap(CVmix_vars, CVmix_kpp_params_user)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for the double diffusion mixing
!  parameterizatiion.
!\\
!\\
!
! !USES:
!  only those used by entire module.

! !INPUT PARAMETERS:
    type(cvmix_kpp_params_type), intent(in), optional, target ::              &
                                           CVmix_kpp_params_user

! !INPUT/OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars

!EOP
!BOC

    call cvmix_put(CVmix_vars, 'kpp_transport', 0.0_cvmix_r8)
    call cvmix_coeffs_kpp(CVmix_vars%diff_iface, CVmix_vars%visc_iface,       &
                          CVmix_vars%zw_iface, CVmix_vars%zt,                 &
                          CVmix_vars%OBL_depth, CVmix_vars%kOBL_depth,        &
                          CVmix_vars%kpp_transport_iface,                     &
                          CVmix_vars%surf_fric, CVmix_vars%surf_buoy,         &
                          CVmix_kpp_params_user)

!EOC

  end subroutine cvmix_coeffs_kpp_wrap

!BOP

! !IROUTINE: cvmix_coeffs_kpp_low
! !INTERFACE:

  subroutine cvmix_coeffs_kpp_low(diff, visc, zw_iface, zt_cntr, OBL_depth,   &
                                  kOBL_depth, nonlocal, surf_fric, surf_buoy, &
                                  CVmix_kpp_params_user)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for the double diffusion mixing
!  parameterizatiion.
!\\
!\\
!
! !USES:
!  only those used by entire module.

! !INPUT PARAMETERS:
    type(cvmix_kpp_params_type), intent(in), optional, target ::              &
                                           CVmix_kpp_params_user
    real(cvmix_r8), dimension(:),   intent(in) :: zw_iface, zt_cntr
    real(cvmix_r8),                 intent(in) :: OBL_depth, surf_fric,       &
                                                  surf_buoy, kOBL_depth

! !INPUT/OUTPUT PARAMETERS:
    real(cvmix_r8), dimension(:,:), intent(inout) :: diff, nonlocal
    real(cvmix_r8), dimension(:),   intent(inout) :: visc

!EOP
!BOC

    ! Local variables
    type(cvmix_kpp_params_type), pointer      :: CVmix_kpp_params_in

    ! OBL_diff and OBL_visc are the diffusivity and viscosity in the whole OBL
    real(cvmix_r8), dimension(:,:), allocatable :: OBL_diff
    real(cvmix_r8), dimension(:),   allocatable :: OBL_visc

    ! diff_ktup and visc_ktup are the enhanced diffusivity and viscosity values
    ! at the deepest cell center above OBL_depth. Rest are intermediary
    ! variables needed to compute diff_ktup and visc_ktup
    real(cvmix_r8), dimension(2) :: diff_ktup
    real(cvmix_r8)               :: visc_ktup
    real(cvmix_r8)               :: sigma_ktup, wm_ktup, ws_ktup

    ! enh_diff and enh_visc are the enhanced diffusivity and viscosity values
    ! at the interface nearest OBL_depth
    real(cvmix_r8), dimension(2) :: enh_diff
    real(cvmix_r8)               :: enh_visc
    real(cvmix_r8)               :: delta, omd

    real(cvmix_r8), dimension(:), allocatable :: sigma, w_m, w_s
    real(cvmix_r8), dimension(4,3)            :: shape_coeffs
    real(cvmix_r8), dimension(3) :: Gat1, DGat1, GatS, visc_at_OBL, dvisc_OBL
    real(cvmix_r8) :: wm_OBL, ws_OBL, second_term

    ! Constants from params
    real(cvmix_r8) :: Cstar, vonkar, c_s, surf_layer_ext
    integer :: interp_type2

    integer :: nlev_p1, nlev, kw, i
    logical :: lstable
    integer :: ktup, & ! kt index of cell center above OBL_depth
               kwup    ! kw index of iface above OBL_depth (= kt index of
                       ! cell containing OBL_depth)

    CVmix_kpp_params_in => CVmix_kpp_params_saved
    if (present(CVmix_kpp_params_user)) then
      CVmix_kpp_params_in => CVmix_kpp_params_user
    end if
    interp_type2   = CVmix_kpp_params_in%interp_type2
    vonkar         = cvmix_get_kpp_real('vonkarman', CVmix_kpp_params_in)
    Cstar          = cvmix_get_kpp_real('Cstar', CVmix_kpp_params_in)
    surf_layer_ext = cvmix_get_kpp_real('surf_layer_ext', CVmix_kpp_params_in)
    c_s            = cvmix_get_kpp_real('c_s', CVmix_kpp_params_in)


    nlev_p1 = size(zw_iface)
    nlev    = size(zt_cntr)
    allocate(sigma(nlev_p1), w_m(nlev_p1), w_s(nlev_p1))
    sigma = -zw_iface/OBL_depth

    kwup = floor(kOBL_depth)
    ktup = nint(kOBL_depth)-1

    ! Allocate OBL_diff and OBL_visc
    allocate(OBL_diff(kwup,2), OBL_visc(kwup))
    OBL_diff = 0.0_cvmix_r8
    OBL_visc = 0.0_cvmix_r8

    ! Stability => positive surface buoyancy flux
    lstable = (surf_buoy.gt.0.0_cvmix_r8)

    ! (1) Compute turbulent velocity scales in column and at OBL_depth. Per
    !     
    call cvmix_kpp_compute_turbulent_scales(sigma, OBL_depth, surf_buoy,      &
                                            surf_fric, w_m, w_s)
    call cvmix_kpp_compute_turbulent_scales(1.0_cvmix_r8, OBL_depth,          &
                                            surf_buoy, surf_fric, wm_OBL,     &
                                            ws_OBL)

    ! (2) Compute G(1) and G'(1) for three cases:
    !     i) temperature diffusivity
    !     ii) other tracers diffusivity
    !     iii) viscosity
    if (kwup.eq.1) then
      visc_at_OBL(1) = compute_nu_at_OBL_depth(interp_type2, (/zw_iface(kwup),&
                         zw_iface(kwup+1)/), (/diff(kwup,1), diff(kwup+1,1)/),&
                         OBL_depth, dnu_dz=dvisc_OBL(1))
      visc_at_OBL(2) = compute_nu_at_OBL_depth(interp_type2, (/zw_iface(kwup),&
                         zw_iface(kwup+1)/), (/diff(kwup,2), diff(kwup+1,2)/),&
                         OBL_depth, dnu_dz=dvisc_OBL(2))
      visc_at_OBL(3) = compute_nu_at_OBL_depth(interp_type2, (/zw_iface(kwup),&
                         zw_iface(kwup+1)/), (/visc(kwup), visc(kwup+1)/),    &
                         OBL_depth, dnu_dz=dvisc_OBL(3))
    else
      visc_at_OBL(1) = compute_nu_at_OBL_depth(interp_type2, (/zw_iface(kwup),&
                         zw_iface(kwup+1)/), (/diff(kwup,1), diff(kwup+1,1)/),&
                         OBL_depth, zw_iface(kwup-1), diff(kwup-1,1),         &
                         dvisc_OBL(1)) 
      visc_at_OBL(2) = compute_nu_at_OBL_depth(interp_type2, (/zw_iface(kwup),&
                         zw_iface(kwup+1)/), (/diff(kwup,2), diff(kwup+1,2)/),&
                         OBL_depth, zw_iface(kwup-1), diff(kwup-1,2),         &
                         dvisc_OBL(2)) 
      visc_at_OBL(3) = compute_nu_at_OBL_depth(interp_type2, (/zw_iface(kwup),&
                         zw_iface(kwup+1)/), (/visc(kwup), visc(kwup+1)/),    &
                         OBL_depth, zw_iface(kwup-1), visc(kwup-1),           &
                         dvisc_OBL(3))
    end if
    Gat1(1) = visc_at_OBL(1)/(OBL_depth*ws_OBL)
    Gat1(2) = visc_at_OBL(2)/(OBL_depth*ws_OBL)
    Gat1(3) = visc_at_OBL(3)/(OBL_depth*wm_OBL)
    if (CVmix_kpp_params_in%lnoDGat1) then
      DGat1   = 0.0_cvmix_r8
    else
      DGat1(1) = -dvisc_OBL(1)/ws_OBL
      DGat1(2) = -dvisc_OBL(2)/ws_OBL
      DGat1(3) = -dvisc_OBL(3)/wm_OBL
      if (lstable) then
        second_term = real(5,cvmix_r8)*surf_buoy/(surf_fric**4)
        DGat1(1) = DGat1(1) + second_term*visc_at_OBL(1)
        DGat1(2) = DGat1(2) + second_term*visc_at_OBL(2)
        DGat1(3) = DGat1(3) + second_term*visc_at_OBL(3)
      end if
    end if

    ! (3) Compute coefficients of shape function
    do i=1,3
      call cvmix_kpp_compute_shape_function_coeffs(Gat1(i), DGat1(i),         &
                                                   shape_coeffs(:,i))
    end do

    ! (4) Compute diffusivities and viscosity in ocean boundary layer and at
    !     the z = zt_cntr(ktup) [z coordinate of last cell center still in the
    !     OBL]. Also compute the non-local transport terms (see note about
    !     what is actually stored in "nonlocal")
    nonlocal = 0.0_cvmix_r8
    sigma_ktup = -zt_cntr(ktup)/OBL_depth
    call cvmix_kpp_compute_turbulent_scales(sigma_ktup, OBL_depth, surf_buoy, &
                                            surf_fric, wm_ktup, ws_ktup)
    do i=1,3
      GatS(i) = cvmix_math_evaluate_cubic(shape_coeffs(:,i), sigma_ktup)
    end do
    diff_ktup(1)   = OBL_depth * ws_ktup * GatS(1)
    diff_ktup(2)   = OBL_depth * ws_ktup * GatS(2)
    visc_ktup      = OBL_depth * wm_ktup * GatS(3)
    do kw=1,kwup
      do i=1,3
        GatS(i) = cvmix_math_evaluate_cubic(shape_coeffs(:,i), sigma(kw))
      end do
      OBL_diff(kw,1)   = OBL_depth * w_s(kw) * GatS(1)
      OBL_diff(kw,2)   = OBL_depth * w_s(kw) * GatS(2)
      OBL_visc(kw)     = OBL_depth * w_m(kw) * GatS(3)
      if (.not.lstable) then
        nonlocal(kw,1) = GatS(1)*(Cstar*vonkar*(vonkar*surf_layer_ext*c_s)**  &
                         (real(1,cvmix_r8)/real(3,cvmix_r8)))
        nonlocal(kw,2) = GatS(2)*(Cstar*vonkar*(vonkar*surf_layer_ext*c_s)**  &
                         (real(1,cvmix_r8)/real(3,cvmix_r8)))
      end if
    end do

    ! (5) Compute enhanced mixing at interface nearest to OBL_depth
    ! OBL_depth is between the centers of cells ktup and ktup+1. The enhanced
    ! mixing described in Appendix D of LMD94 changes the diffusivity values
    ! at the interface between them (interface ktup+1), but the change depends
    ! on whether OBL_depth is in layer ktup (=> kwup = ktup) or layer ktup+1
    ! (=> kwup = ktup+1). I.e. is the interface between them in the OBL or
    ! below it?
    if (ktup.eq.0) then
      ! OBL_depth between surface and first cell center, assume zt_cntr(0)=0
      delta = OBL_depth/(-zt_cntr(ktup+1))
    else
      if (ktup.eq.nlev) then
        ! OBL_depth between bottom cell center and ocean bottom, assume
        ! zt_cntr(ktup+1) = ocn_bottom (which is zw_iface(nlev+1)
        delta = (OBL_depth+zt_cntr(ktup))/(zt_cntr(ktup)-zw_iface(ktup+1))
      else
        delta = (OBL_depth+zt_cntr(ktup))/(zt_cntr(ktup)-zt_cntr(ktup+1))
      end if
    end if
    omd   = 1.0_cvmix_r8 - delta ! omd = one minus delta
    if (ktup.eq.kwup) then
      ! Interface is NOT in the OBL
      ! (a) compute enhanced diffs
      enh_diff(1) = (omd**2)*diff_ktup(1) + (delta**2)*diff(ktup+1,1)
      enh_diff(2) = (omd**2)*diff_ktup(2) + (delta**2)*diff(ktup+1,2)
      enh_visc    = (omd**2)*visc_ktup    + (delta**2)*visc(ktup+1)
      
      ! (b) modify diffusivity (in diff and visc, since OBL_diff and OBL_visc
      !     are not defined at this interface)
      diff(ktup+1,1) = omd*diff(ktup+1,1) + delta*enh_diff(1)
      diff(ktup+1,2) = omd*diff(ktup+1,2) + delta*enh_diff(2)
      visc(ktup+1)   = omd*visc(ktup+1)   + delta*enh_visc
    else
      if (ktup.eq.kwup-1) then
        ! Interface is in the OBL
        ! (a) compute enhanced diffs
        enh_diff(1) = (omd**2)*diff_ktup(1) + (delta**2)*OBL_diff(ktup+1,1)
        enh_diff(2) = (omd**2)*diff_ktup(2) + (delta**2)*OBL_diff(ktup+1,2)
        enh_visc    = (omd**2)*visc_ktup    + (delta**2)*OBL_visc(ktup+1)
      
        ! (b) modify diffusivity (in diff and visc, since OBL_diff and OBL_visc
        !     are not defined at this interface)
        OBL_diff(ktup+1,1) = omd*diff(ktup+1,1) + delta*enh_diff(1)
        OBL_diff(ktup+1,2) = omd*diff(ktup+1,2) + delta*enh_diff(2)
        OBL_visc(ktup+1)   = omd*visc(ktup+1)   + delta*enh_visc
      else
        print*, "ERROR: kwup should be either ktup or ktup+1!"
        deallocate(sigma, w_m, w_s)
        deallocate(OBL_diff, OBL_visc)
        stop 1
      end if
    end if

    ! (6) Combine interior and boundary coefficients
    diff(1:kwup,:) = OBL_diff
    visc(1:kwup) = OBL_visc

    ! Clean up memory
    deallocate(sigma, w_m, w_s)
    deallocate(OBL_diff, OBL_visc)

!EOC
  end subroutine cvmix_coeffs_kpp_low

!BOP

! !IROUTINE: cvmix_put_kpp_real
! !INTERFACE:

  subroutine cvmix_put_kpp_real(varname, val, CVmix_kpp_params_user)

! !DESCRIPTION:
!  Write a real value into a cvmix\_kpp\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    real(cvmix_r8),   intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_kpp_params_type), intent(inout), target, optional ::           &
                                              CVmix_kpp_params_user

!EOP
!BOC

    type(cvmix_kpp_params_type), pointer :: CVmix_kpp_params_out

    CVmix_kpp_params_out => CVmix_kpp_params_saved
    if (present(CVmix_kpp_params_user)) then
      CVmix_kpp_params_out => CVmix_kpp_params_user
    end if

    select case (trim(varname))
      case ('Ri_crit')
        CVmix_kpp_params_out%Ri_crit = val
      case ('vonkarman')
        CVmix_kpp_params_out%vonkarman = val
      case ('Cstar')
        CVmix_kpp_params_out%Cstar = val
      case ('zeta_m')
        CVmix_kpp_params_out%zeta_m = val
      case ('zeta_s')
        CVmix_kpp_params_out%zeta_s = val
      case ('a_m')
        CVmix_kpp_params_out%a_m = val
      case ('a_s')
        CVmix_kpp_params_out%a_s = val
      case ('c_m')
        CVmix_kpp_params_out%c_m = val
      case ('c_s')
        CVmix_kpp_params_out%c_s = val
      case ('surf_layer_ext')
        CVmix_kpp_params_out%surf_layer_ext = val
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
    end select

!EOC

  end subroutine cvmix_put_kpp_real

!BOP

! !IROUTINE: cvmix_put_kpp_int
! !INTERFACE:

  subroutine cvmix_put_kpp_int(CVmix_kpp_params, varname, val)

! !DESCRIPTION:
!  Write an integer value into a cvmix\_kpp\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    integer,          intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_kpp_params_type), intent(inout) :: CVmix_kpp_params
!EOP
!BOC

    select case (trim(varname))
      case ('interp_type')
        CVmix_kpp_params%interp_type = val
      case ('interp_type2')
        CVmix_kpp_params%interp_type2 = val
      case DEFAULT
        call cvmix_put_kpp(varname, real(val, cvmix_r8), CVmix_kpp_params)
    end select

!EOC

  end subroutine cvmix_put_kpp_int

!BOP

! !IROUTINE: cvmix_put_kpp_logical
! !INTERFACE:

  subroutine cvmix_put_kpp_logical(CVmix_kpp_params, varname, val)

! !DESCRIPTION:
!  Write a Boolean value into a cvmix\_kpp\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    logical,          intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_kpp_params_type), intent(inout) :: CVmix_kpp_params
!EOP
!BOC

    select case (trim(varname))
      case ('lEkman')
        CVmix_kpp_params%lEkman = val
      case ('lMonOb')
        CVmix_kpp_params%lMonOb = val
      case ('lnoDGat1')
        CVmix_kpp_params%lnoDGat1 = val
      case ('lavg_N_or_Nsqr')
        CVmix_kpp_params%lavg_N_or_Nsqr = val
      case DEFAULT
        print*, "ERROR: ", trim(varname), " is not a boolean variable!"
        stop 1
    end select

!EOC

  end subroutine cvmix_put_kpp_logical

!BOP

! !IROUTINE: cvmix_get_kpp_real
! !INTERFACE:

  function cvmix_get_kpp_real(varname, CVmix_kpp_params_user)

! !DESCRIPTION:
!  Return the real value of a cvmix\_kpp\_params\_type variable.
!  NOTE: This function is not efficient and is only for infrequent
!  queries of ddiff parameters, such as at initialization.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*),                              intent(in) :: varname
    type(cvmix_kpp_params_type), optional, target, intent(in) ::              &
                                           CVmix_kpp_params_user

! !OUTPUT PARAMETERS:
    real(cvmix_r8) :: cvmix_get_kpp_real
!EOP
!BOC

    type(cvmix_kpp_params_type), pointer :: CVmix_kpp_params_in

    CVmix_kpp_params_in => CVmix_kpp_params_saved
    if (present(CVmix_kpp_params_user)) then
      CVmix_kpp_params_in => CVmix_kpp_params_user
    end if

    cvmix_get_kpp_real = 0.0_cvmix_r8
    select case (trim(varname))
      case ('Ri_crit')
        cvmix_get_kpp_real = CVmix_kpp_params_in%Ri_crit
      case ('vonkarman')
        cvmix_get_kpp_real = CVmix_kpp_params_in%vonkarman
      case ('Cstar')
        cvmix_get_kpp_real = CVmix_kpp_params_in%Cstar
      case ('zeta_m')
        cvmix_get_kpp_real = CVmix_kpp_params_in%zeta_m
      case ('zeta_s')
        cvmix_get_kpp_real = CVmix_kpp_params_in%zeta_s
      case ('a_m')
        cvmix_get_kpp_real = CVmix_kpp_params_in%a_m
      case ('a_s')
        cvmix_get_kpp_real = CVmix_kpp_params_in%a_s
      case ('c_m')
        cvmix_get_kpp_real = CVmix_kpp_params_in%c_m
      case ('c_s')
        cvmix_get_kpp_real = CVmix_kpp_params_in%c_s
      case ('surf_layer_ext')
        cvmix_get_kpp_real = CVmix_kpp_params_in%surf_layer_ext
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
    end select

!EOC

  end function cvmix_get_kpp_real

!BOP

! !IROUTINE: cvmix_kpp_compute_OBL_depth_low
! !INTERFACE:

  subroutine cvmix_kpp_compute_OBL_depth_low(Ri_bulk, zw_iface, OBL_depth,    &
                                             kOBL_depth, zt_cntr, surf_fric,  &
                                             surf_buoy, Coriolis,             &
                                             CVmix_kpp_params_user)

! !DESCRIPTION:
!  Computes the depth of the ocean boundary layer (OBL) for a given column
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    real(cvmix_r8), dimension(:),                   intent(in) :: Ri_bulk
    real(cvmix_r8), dimension(:),           target, intent(in) :: zw_iface,   &
                                                                  zt_cntr
    real(cvmix_r8),               optional,         intent(in) :: surf_fric,  &
                                                                  surf_buoy,  &
                                                                  Coriolis
    type(cvmix_kpp_params_type),  optional, target, intent(in) ::             &
                                            CVmix_kpp_params_user

! !OUTPUT PARAMETERS:
    real(cvmix_r8),               intent(out) :: OBL_depth, kOBL_depth

!EOP
!BOC

    ! Local variables
    real(kind=cvmix_r8), dimension(:), pointer :: depth
    real(kind=cvmix_r8), dimension(4)          :: coeffs
    real(kind=cvmix_r8) :: Ekman, MoninObukhov, OBL_Limit
    integer             :: nlev, k
    logical             :: lstable

    type(cvmix_kpp_params_type), pointer :: CVmix_kpp_params_in

    CVmix_kpp_params_in => CVmix_kpp_params_saved
    if (present(CVmix_kpp_params_user)) then
      CVmix_kpp_params_in => CVmix_kpp_params_user
    end if

    ! Error checks
    ! (1) if using Ekman length, need to pass surf_fric and Coriolis
    if ((.not.(present(surf_fric).and.present(Coriolis))).and.                &
        CVmix_kpp_params_in%lMonOb) then
      print*, "ERROR: must pass surf_fric and Coriolis if you want to ",      &
                "compute Ekman length"
      stop 1
    end if

    ! (2) if using Monin-Obukhov length, need to pass surf_fric and surf_buoy
    if ((.not.(present(surf_fric).and.present(surf_buoy))).and.               &
        CVmix_kpp_params_in%lMonOb) then
      print*, "ERROR: must pass surf_fric and surf_buoy if you want to ",     &
                "compute Monin-Obukhov length"
      stop 1
    end if

    ! (3) Ri_bulk needs to be either the size of zw_iface or zt_cntr
    nlev = size(zw_iface)-1
    if (size(Ri_bulk).eq.nlev) then
      if (size(zt_cntr).ne.nlev) then
        print*, "ERROR: zt_cntr must have length nlev!"
        stop 1
      end if
      depth => zt_cntr
    else
      if (size(Ri_bulk).eq.nlev+1) then
        depth => zw_iface
      else
        print*, "ERROR: Ri_bulk must have size nlev or nlev+1!"
        stop 1
      end if
    end if

    ! if lEkman = .true., OBL_depth must be between the surface and the Ekman
    ! depth. Similarly, if lMonOb = .true., OBL_depth must be between the
    ! surface and the Monin-Obukhov depth
    OBL_limit  = abs(zt_cntr(nlev))

    ! Since depth gets more negative as you go deeper, that translates into
    ! OBL_depth = max(computed depth, Ekman depth, M-O depth)
    ! (MNL: change this when we make OBL_depth positive-down!)
    if (CVmix_kpp_params_in%lEkman) then
      if (Coriolis.eq.0.0_cvmix_r8) then
        ! Rather than divide by zero, set Ekman depth to ocean bottom
        Ekman = abs(zt_cntr(nlev))
      else
        Ekman = 0.7_cvmix_r8*surf_fric/abs(Coriolis)
      end if
      OBL_limit = min(OBL_limit, Ekman)
    end if

    if (CVmix_kpp_params_in%lMonOb) then
      ! Column is stable if surf_buoy > 0
      lstable = (surf_buoy.gt.0.0_cvmix_r8)

      if (lstable) then
        MoninObukhov = surf_fric**3/(surf_buoy*cvmix_get_kpp_real('vonkarman',&
                                                     CVmix_kpp_params_in))
      else
        MoninObukhov = abs(zt_cntr(nlev))
      end if
      OBL_limit = min(OBL_limit, MoninObukhov)
    end if

    ! Interpolation Step
    ! (1) Find k such that Ri_bulk at level k+1 > Ri_crit
    do k=1,size(Ri_bulk)-1
      if (Ri_bulk(k+1).gt.CVmix_kpp_params_in%ri_crit) &
        exit
    end do

    if (k.eq.size(Ri_bulk)) then
      OBL_depth = abs(OBL_limit)
    else
      if (k.eq.1) then
        call cvmix_math_poly_interp(coeffs, CVmix_kpp_params_in%interp_type,  &
                               depth(k:k+1), Ri_bulk(k:k+1))
      else
        call cvmix_math_poly_interp(coeffs, CVmix_kpp_params_in%interp_type,  &
                               depth(k:k+1), Ri_bulk(k:k+1), depth(k-1),      &
                               Ri_bulk(k-1))
      end if
      coeffs(1) = coeffs(1)-CVmix_kpp_params_in%ri_crit

      OBL_depth = -cvmix_math_cubic_root_find(coeffs,                         &
                                         0.5_cvmix_r8*(depth(k)+depth(k+1)))

      ! Note: maybe there are times when we don't need to do the interpolation
      !       because we know OBL_depth will equal OBL_limit?
      OBL_depth = min(OBL_depth, OBL_limit)
    end if

    kOBL_depth = cvmix_kpp_compute_kOBL_depth(zw_iface, zt_cntr, OBL_depth)

!EOC

  end subroutine cvmix_kpp_compute_OBL_depth_low

!BOP

! !IROUTINE: cvmix_kpp_compute_kOBL_depth
! !INTERFACE:

  function cvmix_kpp_compute_kOBL_depth(zw_iface, zt_cntr, OBL_depth)

! !DESCRIPTION:
!  Computes the index of the level and interface above OBL_depth. The index is
!  stored as a real number, and the integer index can be solved for in the
!  following way:\\
!    \verb|kt| = index of cell center above OBL_depth = \verb|nint(kOBL_depth)-1|
!    \verb|kw| = index of interface above OBL_depth = \verb|floor(kOBL_depth)|
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    real(cvmix_r8), dimension(:), intent(in) :: zw_iface, zt_cntr
    real(cvmix_r8), intent(in)  :: OBL_depth

! !OUTPUT PARAMETERS:
    real(cvmix_r8) :: cvmix_kpp_compute_kOBL_depth

!EOP
!BOC

    ! Local variables
    integer :: kw, nlev

    nlev = size(zt_cntr)
    if (size(zw_iface).ne.nlev+1) then
      print*, "ERROR: there should be one more interface z coordinate than ", &
              "cell center coordinate!"
      stop 1
    end if

    ! Initial value = nlev + 0.75 => OBL_depth at center of bottom cell
    cvmix_kpp_compute_kOBL_depth = real(nlev,cvmix_r8)+0.75_cvmix_r8
    do kw=1,nlev
      if (OBL_depth.lt.abs(zw_iface(kw+1))) then
        if (OBL_depth.lt.abs(zt_cntr(kw))) then
          cvmix_kpp_compute_kOBL_depth = real(kw, cvmix_r8)+0.25_cvmix_r8
        else
          cvmix_kpp_compute_kOBL_depth = real(kw, cvmix_r8)+0.75_cvmix_r8
        end if
        exit
      end if
    end do

!EOC

  end function cvmix_kpp_compute_kOBL_depth

!BOP

! !IROUTINE: cvmix_kpp_compute_OBL_depth_wrap
! !INTERFACE:

  subroutine cvmix_kpp_compute_OBL_depth_wrap(CVmix_vars, CVmix_kpp_params_user)

! !DESCRIPTION:
!  Computes the depth of the ocean boundary layer (OBL) for a given column
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    type(cvmix_kpp_params_type), optional, target, intent(in) ::                &
                                           CVmix_kpp_params_user

! !OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars

!EOP
!BOC

    ! Local variables
    real(cvmix_r8) :: lcl_obl_depth, lcl_kobl_depth

    call cvmix_kpp_compute_OBL_depth(CVmix_vars%Rib, CVmix_vars%zw_iface,     &
                                     lcl_obl_depth,  lcl_kobl_depth,          &
                                     CVmix_vars%zt,                           &
                                     CVmix_vars%surf_fric,                    &
                                     CVmix_vars%surf_buoy,                    & 
                                     CVmix_vars%Coriolis,                     &
                                     CVmix_kpp_params_user)
    call cvmix_put(CVmix_vars, 'OBL_depth', lcl_obl_depth)
    call cvmix_put(CVmix_vars, 'kOBL_depth', lcl_kobl_depth)

!EOC

  end subroutine cvmix_kpp_compute_OBL_depth_wrap

!BOP

! !IROUTINE: cvmix_kpp_compute_bulk_Richardson
! !INTERFACE:

  function cvmix_kpp_compute_bulk_Richardson(zt_cntr, delta_buoy_cntr,        &
                                             delta_Vsqr_cntr, Vt_sqr_cntr,    &
                                             ws_cntr, N_iface, Nsqr_iface,    &
                                             CVmix_kpp_params_user)

! !DESCRIPTION:
!  Computes the bulk Richardson number at cell centers. If \verb|Vt_sqr_cntr|
!  is not present, this routine will call \verb|compute_unresolved_shear|,
!  a routine that requires \verb|ws_cntr| and either \verb|N_iface| or
!  \verb|Nsqr_iface|.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    ! * zt_cntr is level-center height (d in LMD94, units: m)
    ! * delta_buoy_cntr is the mean buoyancy estimate over surface layer minus
    !   the level-center buoyancy ( (Br-B(d)) in LMD94, units: m/s^2)
    ! * delta_Vsqr_cntr is the square of the magnitude of the mean velocity
    !   estimate over surface layer minus the level-center velocity
    !   ( |Vr-V(d)|^2 in LMD94, units: m^2/s^2)
    real(cvmix_r8), dimension(:), intent(in) :: zt_cntr, delta_buoy_cntr,     &
                                                delta_Vsqr_cntr
    ! * ws_cntr: w_s (turbulent scale factor) at center of cell (units: m/s)
    ! * N_iface: buoyancy frequency at interfaces (units: 1/s)
    ! * Nsqr_iface: squared buoyancy frequency at interfaces (units: 1/s^2)
    ! * Vt_sqr_cntr: squared unresolved shear term (units m^2/s^2)
    ! See note in description about what values should be passed in
    real(cvmix_r8), dimension(:), intent(in), optional :: ws_cntr, N_iface,   &
                                                          Nsqr_iface,         &
                                                          Vt_sqr_cntr
    type(cvmix_kpp_params_type), intent(in), optional, target ::              &
                                           CVmix_kpp_params_user

! !OUTPUT PARAMETERS:
    real(cvmix_r8), dimension(size(zt_cntr)) ::                               &
                             cvmix_kpp_compute_bulk_Richardson

!EOP
!BOC

    ! Local variables
    ! * unresolved_shear_cntr_sqr is the square of the unresolved level-center
    !   velocity shear (Vt^2(d) in LMD94, units: m^2/s^2)
    real(cvmix_r8), allocatable, dimension(:) :: unresolved_shear_cntr_sqr
    integer        :: kt
    real(cvmix_r8) :: num, denom

    ! Make sure all arguments are same size
    if (any((/size(delta_buoy_cntr), size(delta_Vsqr_cntr)/).ne.              &
        size(zt_cntr))) then
      print*, "ERROR: delta_buoy, delta_vel_sqr, and zt_cntr must all be the",&
              "same size!"
      stop 1
    end if
    allocate(unresolved_shear_cntr_sqr(size(zt_cntr)))
    if (present(Vt_sqr_cntr)) then
      if (size(Vt_sqr_cntr).eq.size(zt_cntr)) then
        unresolved_shear_cntr_sqr = Vt_sqr_cntr
      else
        print*, "ERROR: Vt_sqr_cntr must be the same size as zt_cntr!"
        stop 1
      end if
    else
      if (.not.present(ws_cntr)) then
        print*, "ERROR: you must pass in either Vt_sqr_cntr or ws_cntr!"
        stop 1
      end if
      unresolved_shear_cntr_sqr = cvmix_kpp_compute_unresolved_shear(zt_cntr, &
                                      ws_cntr, N_iface, Nsqr_iface,           &
                                      CVmix_kpp_params_user)
    end if

    do kt=1,size(zt_cntr)
      ! Negative sign because we use positive-up for height
      num   = -zt_cntr(kt)*delta_buoy_cntr(kt)
      denom = delta_Vsqr_cntr(kt) + unresolved_shear_cntr_sqr(kt)
      if (denom.ne.0.0_cvmix_r8) then
        cvmix_kpp_compute_bulk_Richardson(kt) = num/denom
      else
        ! Need a better fudge factor?
        cvmix_kpp_compute_bulk_Richardson(kt) = num*1e10_cvmix_r8
      end if
    end do
    deallocate(unresolved_shear_cntr_sqr)

!EOC

  end function cvmix_kpp_compute_bulk_Richardson

!BOP

! !IROUTINE: cvmix_kpp_compute_turbulent_scales_0d
! !INTERFACE:

  subroutine cvmix_kpp_compute_turbulent_scales_0d(sigma_coord, OBL_depth,    &
                                                   surf_buoy_force,           &
                                                   surf_fric_vel, w_m, w_s,   &
                                                   CVmix_kpp_params_user)

! !DESCRIPTION:
!  Computes the turbulent velocity scales for momentum ($w\_m$) and scalars
!  ($w\_s$) at single $\sigma$ coordinate
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    real(cvmix_r8), intent(in) :: sigma_coord
    real(cvmix_r8), intent(in) :: OBL_depth, surf_buoy_force, surf_fric_vel
    type(cvmix_kpp_params_type), intent(in), optional, target ::              &
                                           CVmix_kpp_params_user

! !OUTPUT PARAMETERS:
    real(cvmix_r8), optional, intent(inout) :: w_m
    real(cvmix_r8), optional, intent(inout) :: w_s

!EOP
!BOC

    ! Local variables
    real(cvmix_r8), dimension(1) :: sigma, lcl_wm, lcl_ws
    logical :: compute_wm, compute_ws

    compute_wm = present(w_m)
    compute_ws = present(w_s)
    sigma(1) = sigma_coord
    if (compute_wm) &
      lcl_wm(1) = w_m
    if (compute_ws) &
      lcl_ws(1) = w_s
    if (compute_wm.and.compute_ws) then
      call cvmix_kpp_compute_turbulent_scales(sigma, OBL_depth,               &
                                              surf_buoy_force, surf_fric_vel, &
                                              w_m = lcl_wm, w_s = lcl_ws,     &
                                  CVmix_kpp_params_user=CVmix_kpp_params_user)
    else
      if (compute_wm) &
        call cvmix_kpp_compute_turbulent_scales(sigma, OBL_depth,             &
                                                surf_buoy_force,surf_fric_vel,&
                                                w_m = lcl_wm,                 &
                                  CVmix_kpp_params_user=CVmix_kpp_params_user)
      if (compute_ws) &
        call cvmix_kpp_compute_turbulent_scales(sigma, OBL_depth,             &
                                                surf_buoy_force,surf_fric_vel,&
                                                w_s = lcl_ws,                 &
                                  CVmix_kpp_params_user=CVmix_kpp_params_user)
    end if

    if (compute_wm) &
      w_m = lcl_wm(1)
    if (compute_ws) &
      w_s = lcl_ws(1)

!EOC

  end subroutine cvmix_kpp_compute_turbulent_scales_0d

!BOP

! !IROUTINE: cvmix_kpp_compute_turbulent_scales_1d
! !INTERFACE:

  subroutine cvmix_kpp_compute_turbulent_scales_1d(sigma_coord, OBL_depth,    &
                                                   surf_buoy_force,           &
                                                   surf_fric_vel, w_m, w_s,   &
                                                   CVmix_kpp_params_user)

! !DESCRIPTION:
!  Computes the turbulent velocity scales for momentum (\verb|w_m|) and scalars
!  (\verb|w_s|) given a 1d array of $\sigma$ coordinates. Note that the
!  turbulent scales are a continuous function, so there is no restriction to
!  only evaluating this routine at interfaces or cell centers. Also, if 
!  $sigma \gt$ \verb|surf_layer_ext| (which is typically 0.1), \verb|w_m| and
!  \verb|w_s| will be evaluated at the latter value.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    real(cvmix_r8), dimension(:), intent(in) :: sigma_coord
    real(cvmix_r8), intent(in) :: OBL_depth, surf_buoy_force, surf_fric_vel
    type(cvmix_kpp_params_type), intent(in), optional, target ::              &
                                           CVmix_kpp_params_user

! !OUTPUT PARAMETERS:
    real(cvmix_r8), optional, dimension(:), intent(inout) :: w_m
    real(cvmix_r8), optional, dimension(:), intent(inout) :: w_s

!EOP
!BOC

    ! Local variables
    integer :: n_sigma, kw
    logical :: compute_wm, compute_ws
    real(cvmix_r8), allocatable, dimension(:) :: zeta
    real(cvmix_r8) :: vonkar, surf_layer_ext
    type(cvmix_kpp_params_type), pointer :: CVmix_kpp_params_in

    n_sigma = size(sigma_coord)

    CVmix_kpp_params_in => CVmix_kpp_params_saved
    if (present(CVmix_kpp_params_user)) then
      CVmix_kpp_params_in => CVmix_kpp_params_user
    end if

    compute_wm = present(w_m)
    compute_ws = present(w_s)
    vonkar = cvmix_get_kpp_real('vonkarman', CVmix_kpp_params_in)
    surf_layer_ext = cvmix_get_kpp_real('surf_layer_ext', CVmix_kpp_params_in)

    if (surf_fric_vel.ne.0.0_cvmix_r8) then
      allocate(zeta(n_sigma))
      do kw=1,n_sigma
        ! compute scales at sigma if sigma < surf_layer_ext, otherwise compute
        ! at surf_layer_ext
        zeta(kw) = min(surf_layer_ext, sigma_coord(kw)) * OBL_depth *         &
                   surf_buoy_force*vonkar/(surf_fric_vel**3)
      end do

      if (compute_wm) then
        if (size(w_m).ne.n_sigma) then
          print*, "ERROR: sigma_coord and w_m must be same size!"
          deallocate(zeta)
          stop 1
        end if
        w_m(1) = compute_phi_inv(zeta(1), CVmix_kpp_params_in, lphi_m=.true.)*&
                 vonkar*surf_fric_vel
        do kw=2,n_sigma
          if (zeta(kw).eq.zeta(kw-1)) then
            w_m(kw) = w_m(kw-1)
          else
            w_m(kw) = vonkar*surf_fric_vel*compute_phi_inv(zeta(kw),          &
                                           CVmix_kpp_params_in, lphi_m=.true.)
          end if
        end do
      end if

      if (compute_ws) then
        if (size(w_s).ne.n_sigma) then
          print*, "ERROR: sigma_coord and w_s must be same size!"
          deallocate(zeta)
          stop 1
        end if
        w_s(1) = compute_phi_inv(zeta(1), CVmix_kpp_params_in, lphi_s=.true.)*&
                 vonkar*surf_fric_vel
        do kw=2,n_sigma
          if (zeta(kw).eq.zeta(kw-1)) then
            w_s(kw) = w_s(kw-1)
          else
            w_s(kw) = vonkar*surf_fric_vel*compute_phi_inv(zeta(kw),          &
                                           CVmix_kpp_params_in, lphi_s=.true.)
          end if
        end do
      end if

      deallocate(zeta)

    else ! surf_fric_vel = 0
      if (compute_wm) then
        if (size(w_m).ne.n_sigma) then
          print*, "ERROR: sigma_coord and w_m must be same size!"
          stop 1
        end if

        if (surf_buoy_force.ge.0.0_cvmix_r8) then
          ! Stable regime with surf_fric_vel = 0 => w_m = 0
          w_m = 0.0_cvmix_r8
        else
          ! Unstable forcing, Eq. (B1c) reduces to following
          do kw=1,n_sigma
            w_m(kw) = cvmix_get_kpp_real('c_m', CVmix_kpp_params_in) *        &
                      min(surf_layer_ext, sigma_coord(kw)) * vonkar *         &
                      surf_buoy_force
            if (w_m(kw).lt.0.0_cvmix_r8) then
              w_m(kw) = (-w_m(kw))**(real(1,cvmix_r8)/real(3,cvmix_r8))*vonkar
            else
              w_m(kw) = -(w_m(kw))**(real(1,cvmix_r8)/real(3,cvmix_r8))*vonkar
            end if
          end do
        end if ! surf_buoy_force >= 0
      end if ! compute_wm

      if (compute_ws) then
        if (size(w_s).ne.n_sigma) then
          print*, "ERROR: sigma_coord and w_s must be same size!"
          stop 1
        end if

        if (surf_buoy_force.ge.0.0_cvmix_r8) then
          ! Stable regime with surf_fric_vel = 0 => w_s = 0
          w_s = 0.0_cvmix_r8
        else
          ! Unstable forcing, Eq. (B1c) reduces to following
          do kw=1,n_sigma
            w_s(kw) = cvmix_get_kpp_real('c_s', CVmix_kpp_params_in) *        &
                      min(surf_layer_ext, sigma_coord(kw)) * vonkar *         &
                      surf_buoy_force
            if (w_s(kw).lt.0) then
              w_s(kw) = (-w_s(kw))**(real(1,cvmix_r8)/real(3,cvmix_r8))*vonkar
            else
              w_s(kw) = -(w_s(kw))**(real(1,cvmix_r8)/real(3,cvmix_r8))*vonkar
            end if
          end do
        end if ! surf_buoy_force >= 0
      end if ! compute_ws
    end if ! surf_fric_vel != 0

!EOC

  end subroutine cvmix_kpp_compute_turbulent_scales_1d

!BOP

! !IROUTINE: cvmix_kpp_compute_unresolved_shear
! !INTERFACE:

  function cvmix_kpp_compute_unresolved_shear(zt_cntr, ws_cntr, N_iface,      &
                                            Nsqr_iface, CVmix_kpp_params_user)

! !DESCRIPTION:
!  Computes the square of the unresolved shear ($V_t^2$ in Eq. (23) of LMD94)
!  at cell centers. Note that you must provide either the buoyancy frequency
!  or its square at cell interfaces, this routine by default will use the
!  lower cell interface value as the cell center, but you can instead take
!  an average of the top and bottom interface values by setting
!  lavg_N_or_Nsqr = .true. in cvmix_kpp_init(). If you pass in Nsqr then
!  negative values are assumed to be zero (default POP behavior)
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    ! zt_cntr: height at center of cell (units: m)
    ! ws_cntr: w_s (turbulent scale factor) at center of cell (units: m/s)
    real(cvmix_r8), dimension(:), intent(in) :: zt_cntr,  ws_cntr
    ! N_iface: buoyancy frequency at cell interfaces (units: 1/s)
    ! Nsqr_iface: squared buoyancy frequency at cell interfaces (units: 1/s^2)
    ! note that you must provide exactly one of these two inputs!
    real(cvmix_r8), dimension(:), intent(in), optional :: N_iface, Nsqr_iface
    type(cvmix_kpp_params_type),  intent(in), optional, target ::             &
                                           CVmix_kpp_params_user

! !OUTPUT PARAMETERS:
    real(cvmix_r8), dimension(size(zt_cntr)) ::                               &
                             cvmix_kpp_compute_unresolved_shear

!EOP
!BOC

    ! Local variables
    integer :: kt, nlev
    real(cvmix_r8) :: Cv, Vtc
    ! N_cntr: buoyancy frequency at cell centers, derived from either N_iface
    !        or Nsqr_iface (units: 1/s)
    real(cvmix_r8), dimension(:), allocatable :: N_cntr
    type(cvmix_kpp_params_type), pointer :: CVmix_kpp_params_in

    nlev = size(zt_cntr)
    if (size(ws_cntr).ne.nlev) then
      print*, "ERROR: zt_cntr and ws_cntr must be same size"
      stop 1
    end if

    if (present(N_iface).and.present(Nsqr_iface)) then
      print*, "ERROR: you must provide N_iface OR Nsqr_iface, can not send",  &
              "both!"
      stop 1
    end if

    if (.not.(present(N_iface).or.present(Nsqr_iface))) then
    end if

    CVmix_kpp_params_in => CVmix_kpp_params_saved
    if (present(CVmix_kpp_params_user)) then
      CVmix_kpp_params_in => CVmix_kpp_params_user
    end if

    if (present(N_iface)) then
      if (size(N_iface).ne.(nlev+1)) then
        print*, "ERROR: N_iface must have one more element than zt_cntr"
        stop 1
      end if
      allocate(N_cntr(nlev))
      do kt=1,nlev
        if (CVmix_kpp_params_in%lavg_N_or_Nsqr) then
          N_cntr(kt) = 0.5_cvmix_r8*(N_iface(kt)+N_iface(kt+1))
        else
          N_cntr(kt) = N_iface(kt+1)
        end if
      end do
    else
      if (present(Nsqr_iface)) then
        if (size(Nsqr_iface).ne.(nlev+1)) then
          print*, "ERROR: Nsqr_iface must have one more element than zt_cntr"
          stop 1
        end if
        allocate(N_cntr(nlev))
        do kt=1,nlev
          if (CVmix_kpp_params_in%lavg_N_or_Nsqr) then
            N_cntr(kt)=sqrt((max(Nsqr_iface(kt),0.0_cvmix_r8) +               &
                             max(Nsqr_iface(kt+1),0.0_cvmix_r8)) *            &
                             0.5_cvmix_r8)
          else
            N_cntr(kt)=sqrt(max(Nsqr_iface(kt+1),0.0_cvmix_r8))
          end if
        end do
      else
        print*, "ERROR: you must provide N_iface OR Nsqr_iface"
        stop 1
      end if
    end if

    ! From LMD 94, Vtc = sqrt(-beta_T/(c_s*eps))/kappa^2
    Vtc = sqrt(0.2_cvmix_r8/(cvmix_get_kpp_real('c_s', CVmix_kpp_params_in) * &
                cvmix_get_kpp_real('surf_layer_ext', CVmix_kpp_params_in))) / &
          (cvmix_get_kpp_real('vonkarman', CVmix_kpp_params_in)**2)
    do kt=1,nlev
      ! Cv computation comes from Danabasoglu et al., 2006
      if (N_cntr(kt).lt.0.002_cvmix_r8) then
        Cv = 2.1_cvmix_r8-200.0_cvmix_r8*N_cntr(kt)
      else
        Cv = 1.7_cvmix_r8
      end if

      cvmix_kpp_compute_unresolved_shear(kt) = -Cv*Vtc*zt_cntr(kt)*           &
                            N_cntr(kt)*ws_cntr(kt)/CVmix_kpp_params_in%Ri_crit
    end do

    deallocate(N_cntr)

!EOC

  end function cvmix_kpp_compute_unresolved_shear

  function compute_phi_inv(zeta, CVmix_kpp_params_in, lphi_m, lphi_s)

    real(cvmix_r8),              intent(in) :: zeta
    type(cvmix_kpp_params_type), intent(in) :: CVmix_kpp_params_in
    logical, optional,           intent(in) :: lphi_m, lphi_s

    real(cvmix_r8) :: compute_phi_inv

    logical :: lm, ls

    ! If not specifying lphi_m or lphi_s, routine will error out, but
    ! initializing result to 0 removes warning about possibly returning an
    ! un-initialized value
    compute_phi_inv = 0.0_cvmix_r8

    if (present(lphi_m)) then
      lm = lphi_m
    else
      lm = .false.
    end if

    if (present(lphi_s)) then
      ls = lphi_s
    else
      ls = .false.
    end if

    if (lm.eqv.ls) then
      print*, "ERROR: must compute phi_m or phi_s, can not compute both!"
      stop 1
    end if

    if (lm) then
      if (zeta.ge.0.0_cvmix_r8) then
        ! Stable region
        compute_phi_inv = one/(one + real(5,cvmix_r8)*zeta)
      else if (zeta.ge.cvmix_get_kpp_real('zeta_m', CVmix_kpp_params_in)) then
        compute_phi_inv = (one - real(16,cvmix_r8)*zeta)**0.25_cvmix_r8
      else
        compute_phi_inv = (cvmix_get_kpp_real('a_m', CVmix_kpp_params_in) -      &
                          cvmix_get_kpp_real('c_m', CVmix_kpp_params_in)*zeta)** &
                          (one/real(3,cvmix_r8))
      end if
    end if

    if (ls) then
      if (zeta.ge.0.0_cvmix_r8) then
        ! Stable region
        compute_phi_inv = one/(one + real(5,cvmix_r8)*zeta)
      else if (zeta.ge.cvmix_get_kpp_real('zeta_s', CVmix_kpp_params_in)) then
        compute_phi_inv = (one - real(16,cvmix_r8)*zeta)**0.5_cvmix_r8
      else
        compute_phi_inv = (cvmix_get_kpp_real('a_s', CVmix_kpp_params_in) -      &
                          cvmix_get_kpp_real('c_s', CVmix_kpp_params_in)*zeta)** &
                          (one/real(3,cvmix_r8))
      end if
    end if

  end function compute_phi_inv

!BOP

! !IROUTINE: cvmix_kpp_compute_shape_function_coeffs
! !INTERFACE:

  subroutine cvmix_kpp_compute_shape_function_coeffs(GAT1, DGAT1, coeffs)

! !DESCRIPTION:
!  Computes the coefficients of the shape function $G(\sigma) = a_0 + a_1\sigma
!  + a_2\sigma^2 + a_3\sigma^3$, where
!  \begin{eqnarray*}
!    a_0 & = & 0 \\
!    a_1 & = & 1 \\
!    a_2 & = &  3G(1) - G'(1) - 2 \\
!    a_3 & = & -2G(1) + G'(1) + 1
!  \end{eqnarray*}
!  Note that $G(1)$ and $G'(1)$ come from Eq. (18) in Large, et al., and
!  this routine returns coeffs(1:4) = $(/a_0, a_1, a_2, a_3/)$
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    real(cvmix_r8), intent(in) :: GAT1  ! G(1)
    real(cvmix_r8), intent(in) :: DGAT1 ! G'(1)

! !OUTPUT PARAMETERS:
    real(cvmix_r8), dimension(4), intent(inout) :: coeffs

!EOP
!BOC

    coeffs(1) =  0.0_cvmix_r8
    coeffs(2) =  1.0_cvmix_r8
    coeffs(3) =  3.0_cvmix_r8*GAT1 - DGAT1 - 2.0_cvmix_r8
    coeffs(4) = -2.0_cvmix_r8*GAT1 + DGAT1 + 1.0_cvmix_r8

!EOC

  end subroutine cvmix_kpp_compute_shape_function_coeffs

  function compute_nu_at_OBL_depth(interp_type2, layer_depth, layer_nu,       &
                                   OBL_depth, depth_2above, nu_2above, dnu_dz)

! !INPUT PARAMETERS:
    integer,                      intent(in) :: interp_type2
    ! layer_depth = (/depth_above_OBL, depth_below_OBL/)
    ! layer_nu    = nu at these points
    real(cvmix_r8), dimension(2), intent(in) :: layer_depth, layer_nu
    real(cvmix_r8),               intent(in) :: OBL_depth
    ! nu at iface above the iface above OBL_depth (Not needed for linear
    ! interpolation or if OBL_depth is in top level
    real(cvmix_r8), optional,     intent(in) :: depth_2above, nu_2above

! !OUTPUT PARAMETERS:
    real(cvmix_r8), optional, intent(out) :: dnu_dz
    real(cvmix_r8)                        :: compute_nu_at_OBL_depth

    ! Local variables
    real(cvmix_r8), dimension(4) :: coeffs
    real(cvmix_r8) :: dnu_dz_above, dnu_dz_below, dnu_dz_local

    if (interp_type2.eq.CVMIX_KPP_INTERP_POP) then
      ! (1) Interpolate derivatives of nu
      if (present(depth_2above).and.present(nu_2above)) then
        dnu_dz_above = (layer_nu(1)-nu_2above)/(layer_depth(1)-depth_2above)
      else
        dnu_dz_above = 0.0_cvmix_r8
      end if
      dnu_dz_below = (layer_nu(2)-layer_nu(1))/(layer_depth(2)-layer_depth(1))
      call cvmix_math_poly_interp(coeffs, CVMIX_MATH_INTERP_LINEAR,           &
                                  layer_depth, (/dnu_dz_above, dnu_dz_below/))
      ! (2) Evaluate at OBL_depth
      dnu_dz_local = cvmix_math_evaluate_cubic(coeffs, -OBL_depth)
      ! (3) Linear interpolant: slope = value computed in (2) and the line goes
      !     through the point (layer_depth(2), layer_nu(2))
      coeffs = 0.0_cvmix_r8
      coeffs(1) = layer_nu(2) - dnu_dz_local*layer_depth(2)
      coeffs(2) = dnu_dz_local
    else
      call cvmix_math_poly_interp(coeffs, interp_type2, layer_depth, layer_nu,&
           depth_2above, nu_2above)
    end if
    compute_nu_at_OBL_depth = cvmix_math_evaluate_cubic(coeffs, -OBL_depth,   &
                                                        dnu_dz)

  end function compute_nu_at_OBL_depth

end module cvmix_kpp
