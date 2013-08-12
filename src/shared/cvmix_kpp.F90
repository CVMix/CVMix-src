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
  public :: cvmix_kpp_compute_turbulent_scales
  ! These are public for testing, may end up private later
  public :: cvmix_kpp_compute_shape_function_coeffs

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
    real(cvmix_r8) :: Ri_crit      ! Critical Richardson number
                               ! (OBL_depth = point where bulk Ri = Ri_crit)
    real(cvmix_r8) :: vonkarman    ! von Karman constant
    ! For velocity scale function, _m => momentum and _s => scalar (tracer)
    real(cvmix_r8) :: zeta_m       ! parameter for computing vel scale func
    real(cvmix_r8) :: zeta_s       ! parameter for computing vel scale func
    real(cvmix_r8) :: a_m          ! parameter for computing vel scale func
    real(cvmix_r8) :: a_s          ! parameter for computing vel scale func
    real(cvmix_r8) :: c_m          ! parameter for computing vel scale func
    real(cvmix_r8) :: c_s          ! parameter for computing vel scale func
    real(cvmix_r8) :: eps          ! small non-negative val (rec 1e-10)
    integer        :: interp_type  ! type of interpolation used to interpolate
                                   ! bulk Richardson number
    integer        :: interp_type2 ! type of interpolation used to interpolate
                                   ! diff and visc at OBL_depth
    logical        :: lEkman       ! True => compute Ekman depth limit
    logical        :: lMonOb       ! True => compute Monin-Obukhov limit
    logical        :: lnoDGat1     ! True => G'(1) = 0 (shape function)
                                   ! False => compute G'(1) as in LMD94
  end type cvmix_kpp_params_type

!EOP

type(cvmix_kpp_params_type), target :: CVmix_kpp_params_saved

contains

!BOP

! !IROUTINE: cvmix_init_kpp
! !INTERFACE:

  subroutine cvmix_init_kpp(ri_crit, vonkarman, zeta_m, zeta_s, a_m, a_s,     &
                            c_m, c_s, eps, interp_type, interp_type2, lEkman, &
                            lMonOb, lnoDGat1, CVmix_kpp_params_user)

! !DESCRIPTION:
!  Initialization routine for KPP mixing.
!
! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    real(cvmix_r8),   optional :: ri_crit, vonkarman, zeta_m, zeta_s, a_m, &
                                  a_s, c_m, c_s, eps
    character(len=*), optional :: interp_type, interp_type2
    logical,          optional :: lEkman, lMonOb, lnoDGat1

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

    if (present(eps)) then
      call cvmix_put_kpp('eps', eps, CVmix_kpp_params_user)
    else
      call cvmix_put_kpp('eps', 1e-10_cvmix_r8, CVmix_kpp_params_user)
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

    call cvmix_coeffs_kpp(CVmix_vars%diff_iface, CVmix_vars%visc_iface,       &
                          CVmix_vars%zw_iface, CVmix_vars%OBL_depth,          &
                          CVmix_VARS%kOBL_depth, CVmix_vars%surf_fric,        &
                          CVmix_vars%surf_buoy,         &
                          CVmix_kpp_params_user)

!EOC

  end subroutine cvmix_coeffs_kpp_wrap

!BOP

! !IROUTINE: cvmix_coeffs_kpp_low
! !INTERFACE:

  subroutine cvmix_coeffs_kpp_low(diff, visc, zw_iface, OBL_depth, kup,       &
                                  surf_fric, surf_buoy, CVmix_kpp_params_user)

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
    real(cvmix_r8), dimension(:),   intent(in) :: zw_iface
    real(cvmix_r8),                 intent(in)    :: OBL_depth, surf_fric,    &
                                                     surf_buoy
    integer,                        intent(in)    :: kup ! kw index of iface
                                                         ! above OBL_depth

! !INPUT/OUTPUT PARAMETERS:
    real(cvmix_r8), dimension(:,:), intent(inout) :: diff
    real(cvmix_r8), dimension(:),   intent(inout) :: visc

!EOP
!BOC

    ! Local variables
    type(cvmix_kpp_params_type), pointer      :: CVmix_kpp_params_in
    real(cvmix_r8), dimension(:), allocatable :: sigma, w_m, w_s
    real(cvmix_r8), dimension(4,3)            :: shape_coeffs
    real(cvmix_r8), dimension(3) :: Gat1, DGat1, visc_at_OBL, dvisc_OBL
    real(cvmix_r8)               :: wm_OBL, ws_OBL, second_term
    integer :: nlev_p1, kw, i, interp_type2
    logical :: lstable

    CVmix_kpp_params_in => CVmix_kpp_params_saved
    if (present(CVmix_kpp_params_user)) then
      CVmix_kpp_params_in => CVmix_kpp_params_user
    end if
    interp_type2 = CVmix_kpp_params_in%interp_type2

    nlev_p1 = size(visc)
    allocate(sigma(nlev_p1), w_m(nlev_p1), w_s(nlev_p1))
    sigma = zw_iface/OBL_depth

    ! Stability => positive surface buoyancy flux
    lstable = (surf_buoy.gt.0.0_cvmix_r8)

    ! (1) Compute turbulent velocity scales in column and at OBL_depth
    call cvmix_kpp_compute_turbulent_scales(sigma, OBL_depth, surf_buoy,      &
                                            surf_fric, w_m, w_s)
    call cvmix_kpp_compute_turbulent_scales(1.0_cvmix_r8, OBL_depth,          &
                                            surf_buoy, surf_fric, wm_OBL,     &
                                            ws_OBL)

    ! (2) Compute G(1) and G'(1) for three cases:
    !     i) temperature diffusivity
    !     ii) other tracers diffusivity
    !     iii) viscosity
    if (kup.eq.1) then
      visc_at_OBL(1) = compute_nu_at_OBL_depth(interp_type2, (/zw_iface(kup), &
                         zw_iface(kup+1)/), (/diff(kup,1), diff(kup+1,1)/),   &
                         OBL_depth, dnu_dz=dvisc_OBL(1))
      visc_at_OBL(2) = compute_nu_at_OBL_depth(interp_type2, (/zw_iface(kup), &
                         zw_iface(kup+1)/), (/diff(kup,2), diff(kup+1,2)/),   &
                         OBL_depth, dnu_dz=dvisc_OBL(2))
      visc_at_OBL(3) = compute_nu_at_OBL_depth(interp_type2, (/zw_iface(kup), &
                         zw_iface(kup+1)/), (/visc(kup), visc(kup+1)/),       &
                         OBL_depth, dnu_dz=dvisc_OBL(3))
    else
      visc_at_OBL(1) = compute_nu_at_OBL_depth(interp_type2, (/zw_iface(kup), &
                         zw_iface(kup+1)/), (/diff(kup,1), diff(kup+1,1)/),   &
                         OBL_depth, zw_iface(kup+2), diff(kup+2,1),           &
                         dvisc_OBL(1)) 
      visc_at_OBL(2) = compute_nu_at_OBL_depth(interp_type2, (/zw_iface(kup), &
                         zw_iface(kup+1)/), (/diff(kup,2), diff(kup+1,2)/),   &
                         OBL_depth, zw_iface(kup+2), diff(kup+2,2),           &
                         dvisc_OBL(2)) 
      visc_at_OBL(3) = compute_nu_at_OBL_depth(interp_type2, (/zw_iface(kup), &
                         zw_iface(kup+1)/), (/visc(kup), visc(kup+1)/),       &
                         OBL_depth, zw_iface(kup+2), visc(kup+2), dvisc_OBL(3))
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

    ! (4) Compute diffusivities and viscosity in ocean boundary layer
    do kw=1,kup
      diff(kw,1) = -OBL_depth * w_s(kw) *                                      &
                   cvmix_math_evaluate_cubic(shape_coeffs(:,1), sigma(kw))
      diff(kw,2) = -OBL_depth * w_s(kw) *                                      &
                   cvmix_math_evaluate_cubic(shape_coeffs(:,2), sigma(kw))
      visc(kw)   = -OBL_depth * w_m(kw) *                                      &
                   cvmix_math_evaluate_cubic(shape_coeffs(:,3), sigma(kw))
    end do

    ! (5) Compute non-local transport term

    ! (6) Compute enhanced mixing

    ! (7) Combine interior and boundary coefficients + non-local term

    ! Clean up memory
    deallocate(sigma, w_m, w_s)

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
      case ('eps')
        CVmix_kpp_params_out%eps = val
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
      case ('eps')
        cvmix_get_kpp_real = CVmix_kpp_params_in%eps
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
    real(cvmix_r8), dimension(:),           target, intent(in) :: zw_iface
    real(cvmix_r8), dimension(:), optional, target, intent(in) :: zt_cntr
    real(cvmix_r8),               optional,         intent(in) :: surf_fric,  &
                                                                  surf_buoy,  &
                                                                  Coriolis
    type(cvmix_kpp_params_type),  optional, target, intent(in) ::             &
                                            CVmix_kpp_params_user

! !OUTPUT PARAMETERS:
    real(cvmix_r8), intent(out) :: OBL_depth
    integer,        intent(out) :: kOBL_depth

!EOP
!BOC

    ! Local variables
    real(kind=cvmix_r8), dimension(:), pointer :: depth
    real(kind=cvmix_r8), dimension(4)          :: coeffs
    real(kind=cvmix_r8) :: Ekman, MoninObukhov, OBL_Limit
    integer             :: nlev, k, kw
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
      if (.not.present(zt_cntr)) then
        print*, "ERROR: Ri_bulk has length nlev so you must pass zt_cntr"
        stop 1
      end if
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
    OBL_limit  = depth(nlev)

    ! Since depth gets more negative as you go deeper, that translates into
    ! OBL_depth = max(computed depth, Ekman depth, M-O depth)
    ! (MNL: change this when we make OBL_depth positive-down!)
    if (CVmix_kpp_params_in%lEkman) then
      if (Coriolis.eq.0.0_cvmix_r8) then
        ! Rather than divide by zero, set Ekman depth to ocean bottom
        Ekman = depth(nlev)
      else
        Ekman = 0.7_cvmix_r8*surf_fric/Coriolis
      end if
      OBL_limit = max(OBL_limit, Ekman)
    end if

    if (CVmix_kpp_params_in%lMonOb) then
      ! Column is stable if surf_buoy > 0
      lstable = (surf_buoy.gt.0.0_cvmix_r8)

      if (lstable) then
        MoninObukhov = surf_fric**3/(surf_buoy*cvmix_get_kpp_real('vonkarman',&
                                                     CVmix_kpp_params_in))
      else
        MoninObukhov = depth(nlev)
      end if
      OBL_limit = max(OBL_limit, MoninObukhov)
    end if

    ! Interpolation Step
    ! (1) Find k such that Ri_bulk at level k+1 > Ri_crit
    do k=1,size(Ri_bulk)-1
      if (Ri_bulk(k+1).gt.CVmix_kpp_params_in%ri_crit) &
        exit
    end do
    kOBL_depth = k

    if (k.eq.size(Ri_bulk)) then
      OBL_depth = OBL_limit
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

      OBL_depth = cvmix_math_cubic_root_find(coeffs,                          &
                                         0.5_cvmix_r8*(depth(k)+depth(k+1)))

      ! Note: maybe there are times when we don't need to do the interpolation
      !       because we know OBL_depth will equal OBL_limit?
      OBL_depth = max(OBL_depth, OBL_limit)
    end if

    do kw=1,nlev
      if (OBL_depth.gt.zw_iface(kw+1)) then
        kOBL_depth = kw
        exit
      end if
    end do
!EOC

  end subroutine cvmix_kpp_compute_OBL_depth_low

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
    real(cvmix_r8) :: lcl_obl_depth
    integer        :: lcl_kobl_depth

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
!  Computes the turbulent velocity scales for momentum ($w\_m$) and scalars
!  ($w\_s$) given a 1d array of $\sigma$ coordinates
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
    integer :: nlev_p1, kw
    logical :: compute_wm, compute_ws
    real(cvmix_r8), allocatable, dimension(:) :: zeta, zeta_h
    real(cvmix_r8) :: vonkar
    type(cvmix_kpp_params_type), pointer :: CVmix_kpp_params_in

    nlev_p1 = size(sigma_coord)
    allocate(zeta(nlev_p1), zeta_h(nlev_p1))

    CVmix_kpp_params_in => CVmix_kpp_params_saved
    if (present(CVmix_kpp_params_user)) then
      CVmix_kpp_params_in => CVmix_kpp_params_user
    end if

    compute_wm = present(w_m)
    compute_ws = present(w_s)
    vonkar = cvmix_get_kpp_real('vonkarman', CVmix_kpp_params_in)

    zeta_h = sigma_coord*OBL_depth*surf_buoy_force*vonkar

    zeta = zeta_h/(surf_fric_vel**3 +                                         &
           cvmix_get_kpp_real('eps', CVmix_kpp_params_in))

    if (compute_wm) then
      if (size(w_m).ne.nlev_p1) then
        print*, "ERROR: sigma_coord and w_m must be same size!"
        deallocate(zeta, zeta_h)
        stop 1
      end if
      do kw=1,nlev_p1
        if (zeta(kw).ge.0) then
          ! Stable region
          w_m(kw) = vonkar*surf_fric_vel/(real(1,cvmix_r8) + real(5,cvmix_r8)*&
                    zeta(kw))
        else if (zeta(kw).ge.                                                 &
                 cvmix_get_kpp_real('zeta_m', CVmix_kpp_params_in)) then
          w_m(kw) = vonkar*surf_fric_vel*                                     &
                (real(1,cvmix_r8) - real(16,cvmix_r8)*zeta(kw))**0.25_cvmix_r8
        else
          w_m(kw) = vonkar*(cvmix_get_kpp_real('a_m', CVmix_kpp_params_in)*   &
            (surf_fric_vel**3)-cvmix_get_kpp_real('c_m', CVmix_kpp_params_in)*&
            zeta_h(kw))**(real(1,cvmix_r8)/real(3,cvmix_r8))
        end if
      end do
    end if

    if (compute_ws) then
      if (size(w_s).ne.nlev_p1) then
        print*, "ERROR: sigma_coord and w_s must be same size!"
        deallocate(zeta, zeta_h)
        stop 1
      end if
      do kw=1,nlev_p1
        if (zeta(kw).ge.0) then
          ! Stable region
          w_s(kw) = vonkar*surf_fric_vel/(real(1,cvmix_r8) + real(5,cvmix_r8)*&
                    zeta(kw))
        else if (zeta(kw).ge.                                                 &
                 cvmix_get_kpp_real('zeta_s', CVmix_kpp_params_in)) then
          w_s(kw) = vonkar*surf_fric_vel*                                         &
                sqrt(real(1,cvmix_r8) - real(16,cvmix_r8)*zeta(kw))
        else
          w_s(kw) = vonkar*(cvmix_get_kpp_real('a_s', CVmix_kpp_params_in)*       &
            (surf_fric_vel**3)-cvmix_get_kpp_real('c_s', CVmix_kpp_params_in)*&
            zeta_h(kw))**(real(1,cvmix_r8)/real(3,cvmix_r8))
        end if
      end do
    end if

    deallocate(zeta, zeta_h)

!EOC

  end subroutine cvmix_kpp_compute_turbulent_scales_1d

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
      dnu_dz_local = cvmix_math_evaluate_cubic(coeffs, OBL_depth)
      ! (3) Linear interpolant: slope = value computed in (2) and the line goes
      !     through the point (layer_depth(2), layer_nu(2))
      coeffs = 0.0_cvmix_r8
      coeffs(1) = layer_nu(2) - dnu_dz_local*layer_depth(2)
      coeffs(2) = dnu_dz_local
    else
      call cvmix_math_poly_interp(coeffs, interp_type2, layer_depth, layer_nu,&
           depth_2above, nu_2above)
    end if
    compute_nu_at_OBL_depth = cvmix_math_evaluate_cubic(coeffs, OBL_depth,    &
                                                        dnu_dz)

  end function compute_nu_at_OBL_depth

end module cvmix_kpp
