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

  use cvmix_kinds_and_types, only : cvmix_r8,                &
                                    cvmix_data_type
  use cvmix_put_get, only :         cvmix_put

!EOP

  implicit none
  private
  save

!BOP

! !DEFINED PARAMETERS:
  integer, parameter        :: CVMIX_KPP_INTERP_LINEAR      = 1
  integer, parameter        :: CVMIX_KPP_INTERP_QUAD        = 2
  integer, parameter        :: CVMIX_KPP_INTERP_CUBE_SPLINE = 3
  integer, parameter        :: CVMIX_KPP_MAX_NEWTON_ITERS   = 100
  real(cvmix_r8), parameter :: CVMIX_KPP_NEWTON_TOL         = 1.0e-12_cvmix_r8

! !PUBLIC MEMBER FUNCTIONS:

  public :: cvmix_init_kpp
  public :: cvmix_coeffs_kpp
  public :: cvmix_put_kpp
  public :: cvmix_get_kpp_real
  ! These are public for testing, may end up private later
  public :: cvmix_kpp_compute_OBL_depth
  public :: cvmix_kpp_compute_turbulent_scales
  public :: cvmix_kpp_compute_shape_function_coeffs

  interface cvmix_put_kpp
    module procedure cvmix_put_kpp_int
    module procedure cvmix_put_kpp_real
  end interface cvmix_put_kpp

  interface cvmix_kpp_compute_OBL_depth
    module procedure cvmix_kpp_compute_OBL_depth_low
    module procedure cvmix_kpp_compute_OBL_depth_wrap
  end interface cvmix_kpp_compute_OBL_depth

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
    integer        :: interp_type  ! type of iterpolation to use
  end type cvmix_kpp_params_type

!EOP

type(cvmix_kpp_params_type), target :: CVmix_kpp_params_saved

contains

!BOP

! !IROUTINE: cvmix_init_kpp
! !INTERFACE:

  subroutine cvmix_init_kpp(ri_crit, vonkarman, zeta_m, zeta_s, a_m, a_s,     &
                            c_m, c_s, eps, interp_type, CVmix_kpp_params_user)

! !DESCRIPTION:
!  Initialization routine for KPP mixing.
!
! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    real(cvmix_r8),   optional :: ri_crit, vonkarman, zeta_m, zeta_s, a_m, &
                                  a_s, c_m, c_s, eps
    character(len=*), optional :: interp_type

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
      call cvmix_put_kpp(CVmix_kpp_params_out, 'Ri_crit', ri_crit)
    else
      call cvmix_put_kpp(CVmix_kpp_params_out, 'Ri_crit', 0.3_cvmix_r8)
    end if

    if (present(vonkarman)) then
      call cvmix_put_kpp(CVmix_kpp_params_out, 'vonkarman', vonkarman)
    else
      call cvmix_put_kpp(CVmix_kpp_params_out, 'vonkarman', 0.41_cvmix_r8)
    end if

    if (present(zeta_m)) then
      call cvmix_put_kpp(CVmix_kpp_params_out, 'zeta_m', zeta_m)
    else
      call cvmix_put_kpp(CVmix_kpp_params_out, 'zeta_m', -0.2_cvmix_r8)
    end if

    if (present(zeta_s)) then
      call cvmix_put_kpp(CVmix_kpp_params_out, 'zeta_s', zeta_s)
    else
      call cvmix_put_kpp(CVmix_kpp_params_out, 'zeta_s', -1.0_cvmix_r8)
    end if

    if (present(a_m)) then
      call cvmix_put_kpp(CVmix_kpp_params_out, 'a_m', a_m)
    else
      call cvmix_put_kpp(CVmix_kpp_params_out, 'a_m', 1.26_cvmix_r8)
    end if

    if (present(a_s)) then
      call cvmix_put_kpp(CVmix_kpp_params_out, 'a_s', a_s)
    else
      call cvmix_put_kpp(CVmix_kpp_params_out, 'a_s', -28.86_cvmix_r8)
    end if

    if (present(c_m)) then
      call cvmix_put_kpp(CVmix_kpp_params_out, 'c_m', c_m)
    else
      call cvmix_put_kpp(CVmix_kpp_params_out, 'c_m', 8.38_cvmix_r8)
    end if

    if (present(c_s)) then
      call cvmix_put_kpp(CVmix_kpp_params_out, 'c_s', c_s)
    else
      call cvmix_put_kpp(CVmix_kpp_params_out, 'c_s', 98.96_cvmix_r8)
    end if

    if (present(eps)) then
      call cvmix_put_kpp(CVmix_kpp_params_out, 'eps', eps)
    else
      call cvmix_put_kpp(CVmix_kpp_params_out, 'eps', 1e-10_cvmix_r8)
    end if

    if (present(interp_type)) then
      select case (trim(interp_type))
        case ('line', 'linear')
          call cvmix_put_kpp(CVmix_kpp_params_out, 'interp_type', &
                             CVMIX_KPP_INTERP_LINEAR)
        case ('quad', 'quadratic')
          call cvmix_put_kpp(CVmix_kpp_params_out, 'interp_type', &
                             CVMIX_KPP_INTERP_QUAD)
        case ('cube', 'cubic', 'cubic_spline', 'cubic spline')
          call cvmix_put_kpp(CVmix_kpp_params_out, 'interp_type', &
                             CVMIX_KPP_INTERP_CUBE_SPLINE)
        case DEFAULT
          print*, "ERROR: ", trim(interp_type), " is not a valid type of ", &
                  "interpolation!"
          stop 1
      end select
    else
      call cvmix_put_kpp(CVmix_kpp_params_out, 'interp_type', &
                         CVMIX_KPP_INTERP_QUAD)
    end if
!EOC

  end subroutine cvmix_init_kpp

!***********************************************************************
!BOP
! !IROUTINE: cvmix_coeffs_kpp
! !INTERFACE:

  subroutine cvmix_coeffs_kpp(CVmix_vars, CVmix_kpp_params_user)

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

    type(cvmix_kpp_params_type), pointer :: CVmix_kpp_params_in

    CVmix_kpp_params_in => CVmix_kpp_params_saved
    if (present(CVmix_kpp_params_user)) then
      CVmix_kpp_params_in => CVmix_kpp_params_user
    end if

    call cvmix_kpp_compute_OBL_depth(CVmix_vars, CVmix_kpp_params_in)
    CVmix_vars%visc_iface = cvmix_get_kpp_real('Ri_crit', CVmix_kpp_params_in)
    CVmix_vars%diff_iface = cvmix_get_kpp_real('Ri_crit', CVmix_kpp_params_in)

!EOC
  end subroutine cvmix_coeffs_kpp

!BOP

! !IROUTINE: cvmix_put_kpp_real
! !INTERFACE:

  subroutine cvmix_put_kpp_real(CVmix_kpp_params, varname, val)

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
    type(cvmix_kpp_params_type), intent(inout) :: CVmix_kpp_params
!EOP
!BOC

    select case (trim(varname))
      case ('Ri_crit')
        CVmix_kpp_params%Ri_crit = val
      case ('vonkarman')
        CVmix_kpp_params%vonkarman = val
      case ('zeta_m')
        CVmix_kpp_params%zeta_m = val
      case ('zeta_s')
        CVmix_kpp_params%zeta_s = val
      case ('a_m')
        CVmix_kpp_params%a_m = val
      case ('a_s')
        CVmix_kpp_params%a_s = val
      case ('c_m')
        CVmix_kpp_params%c_m = val
      case ('c_s')
        CVmix_kpp_params%c_s = val
      case ('eps')
        CVmix_kpp_params%eps = val
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
      case DEFAULT
        call cvmix_put_kpp(CVmix_kpp_params, varname, real(val, cvmix_r8))
    end select

!EOC

  end subroutine cvmix_put_kpp_int

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

  subroutine cvmix_kpp_compute_OBL_depth_low(Ri_bulk, depth, OBL_depth, &
                                             CVmix_kpp_params_user)

! !DESCRIPTION:
!  Computes the depth of the ocean boundary layer (OBL) for a given column
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    type(cvmix_kpp_params_type), optional, target, intent(in) ::              &
                                           CVmix_kpp_params_user
    real(cvmix_r8), dimension(:), intent(in) :: Ri_bulk, depth

! !OUTPUT PARAMETERS:
    real(cvmix_r8), intent(out) :: OBL_depth

!EOP
!BOC

    ! Local variables
    integer :: nlev, kt, k
    real(kind=cvmix_r8) :: a, b, c, d, det
    real(kind=cvmix_r8), dimension(:,:), allocatable :: Minv
    real(kind=cvmix_r8), dimension(:),   allocatable :: rhs

    type(cvmix_kpp_params_type), pointer :: CVmix_kpp_params_in

    CVmix_kpp_params_in => CVmix_kpp_params_saved
    if (present(CVmix_kpp_params_user)) then
      CVmix_kpp_params_in => CVmix_kpp_params_user
    end if

    nlev = size(Ri_bulk)
    if (nlev.ne.size(depth)) then
      print*, "ERROR: Ri_bulk and depth must be same size!"
      stop 1
    end if

    ! Interpolation Step
    ! (1) Find kt such that Ri_bulk at level kt+1 > Ri_crit
    do kt=1,nlev-1
      if (Ri_bulk(kt+1).ge.CVmix_kpp_params_in%ri_crit) &
        exit
    end do
    if (kt.eq.nlev) then
      print*, "ERROR: Entire column is above the boundary layer!"
      stop 1
    end if

    ! All interpolation assumes form of
    ! y = ax^3 + bx^2 + cx + d
    ! linear => a = b = 0
    ! quad   => a = 0
    a = 0.0_cvmix_r8
    b = 0.0_cvmix_r8
    c = 0.0_cvmix_r8
    d = 0.0_cvmix_r8
    select case (CVmix_kpp_params_in%interp_type)
      case (CVMIX_KPP_INTERP_LINEAR)
        ! Match values at levels kt and kt+1
        print*, "Linear interpolation"
        c = (Ri_bulk(kt+1)-Ri_bulk(kt))/(depth(kt+1)-depth(kt))
        d = Ri_bulk(kt)-c*depth(kt)
      case (CVMIX_KPP_INTERP_QUAD)
        ! Match slope and value at level kt, value at level kt+1
        print*, "Quadratic interpolation"
        ! [ x1^2 x1 1 ][ b ]   [    y1 ]
        ! [ x0^2 x0 1 ][ c ] = [    y0 ]
        ! [  2x0  1 0 ][ d ]   [ slope ]
        !      ^^^
        !       M
        det = -((depth(kt+1)-depth(kt))**2)
        allocate(Minv(3,3))
        allocate(rhs(3))
        rhs(1) = Ri_bulk(kt+1)
        rhs(2) = Ri_bulk(kt)
        if (kt.gt.1) then
          rhs(3) = (Ri_bulk(kt)-Ri_bulk(kt-1))/(depth(kt)-depth(kt-1))
        else
          rhs(3) = 0.0_cvmix_r8
        end if

        Minv(1,1) = -real(1, cvmix_r8)/det
        Minv(1,2) = real(1, cvmix_r8)/det
        Minv(1,3) = -real(1, cvmix_r8)/(depth(kt+1)-depth(kt))
        Minv(2,1) = real(2, cvmix_r8)*depth(kt)/det
        Minv(2,2) = -real(2, cvmix_r8)*depth(kt)/det
        Minv(2,3) = (depth(kt+1)+depth(kt))/(depth(kt+1)-depth(kt))
        Minv(3,1) = -(depth(kt)**2)/det
        Minv(3,2) = depth(kt+1)*(real(2, cvmix_r8)*depth(kt)-depth(kt+1))/det
        Minv(3,3) = -depth(kt+1)*depth(kt)/(depth(kt+1)-depth(kt))

        do k=1,3
          b = b+Minv(1,k)*rhs(k)
          c = c+Minv(2,k)*rhs(k)
          d = d+Minv(3,k)*rhs(k)
        end do
        deallocate(rhs)
        deallocate(Minv)
      case (CVMIX_KPP_INTERP_CUBE_SPLINE)
        ! [ x1^3 x1^2 x1 1 ][ a ]   [     y1 ]
        ! [ x0^3 x0^2 x0 1 ][ b ] = [     y0 ]
        ! [  3x0  2x0  1 0 ][ c ]   [ slope0 ]
        ! [  3x1  2x1  1 0 ][ d ]   [ slope1 ]
        !      ^^^
        !       M
        det = -((depth(kt+1)-depth(kt))**3)
        allocate(Minv(4,4))
        allocate(rhs(4))
        rhs(1) = Ri_bulk(kt+1)
        rhs(2) = Ri_bulk(kt)
        if (kt.gt.1) then
          rhs(3) = (Ri_bulk(kt)-Ri_bulk(kt-1))/(depth(kt)-depth(kt-1))
        else
          rhs(3) = 0.0_cvmix_r8
        end if
        rhs(4) = (Ri_bulk(kt+1)-Ri_bulk(kt))/(depth(kt+1)-depth(kt))

        Minv(1,1) = real(2, cvmix_r8)/det
        Minv(1,2) = -real(2, cvmix_r8)/det
        Minv(1,3) = (depth(kt)-depth(kt+1))/det
        Minv(1,4) = (depth(kt)-depth(kt+1))/det
        Minv(2,1) = -real(3, cvmix_r8)*(depth(kt+1)+depth(kt))/det
        Minv(2,2) = real(3, cvmix_r8)*(depth(kt+1)+depth(kt))/det
        Minv(2,3) = (depth(kt+1)-depth(kt))*(real(2, cvmix_r8)*depth(kt+1)+depth(kt))/det
        Minv(2,4) = (depth(kt+1)-depth(kt))*(real(2, cvmix_r8)*depth(kt)+depth(kt+1))/det
        Minv(3,1) = real(6, cvmix_r8)*depth(kt+1)*depth(kt)/det
        Minv(3,2) = -real(6, cvmix_r8)*depth(kt+1)*depth(kt)/det
        Minv(3,3) = -depth(kt+1)*(depth(kt+1)-depth(kt))*(real(2, cvmix_r8)*depth(kt)+depth(kt+1))/det
        Minv(3,4) = -depth(kt)*(depth(kt+1)-depth(kt))*(real(2, cvmix_r8)*depth(kt+1)+depth(kt))/det
        Minv(4,1) = -(depth(kt)**2)*(real(3, cvmix_r8)*depth(kt+1)-depth(kt))/det
        Minv(4,2) = -(depth(kt+1)**2)*(-real(3, cvmix_r8)*depth(kt)+depth(kt+1))/det
        Minv(4,3) = depth(kt)*(depth(kt+1)**2)*(depth(kt+1)-depth(kt))/det
        Minv(4,4) = depth(kt+1)*(depth(kt)**2)*(depth(kt+1)-depth(kt))/det

        do k=1,4
          a = a+Minv(1,k)*rhs(k)
          b = b+Minv(2,k)*rhs(k)
          c = c+Minv(3,k)*rhs(k)
          d = d+Minv(4,k)*rhs(k)
        end do
        deallocate(rhs)
        deallocate(Minv)
        ! Match slopes and values at levels kt and kt+1
        print*, "Cubic spline interpolation"
    end select
    print*, kt, nlev
    OBL_depth = cubic_root_find((/a,b,c,d-CVmix_kpp_params_in%ri_crit/), &
                                0.5_cvmix_r8*(depth(kt)+depth(kt+1)))

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

    call cvmix_kpp_compute_OBL_depth(CVmix_vars%Rib, CVmix_vars%zt,           &
                                     lcl_obl_depth, CVmix_kpp_params_user)
    call cvmix_put(CVmix_vars, 'OBL_depth', lcl_obl_depth)

!EOC

  end subroutine cvmix_kpp_compute_OBL_depth_wrap

!BOP

! !IROUTINE: cvmix_kpp_compute_turbulent_scales
! !INTERFACE:

  subroutine cvmix_kpp_compute_turbulent_scales(sigma_coord, OBL_depth,       &
                                                surf_buoy_force,              &
                                                surf_fric_vel, w_m, w_s,      &
                                                CVmix_kpp_params_user)

! !DESCRIPTION:
!  Computes the turbulent velocity scales for momentum ($w\_m$) and scalars
!  ($w\_s$)
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

  end subroutine cvmix_kpp_compute_turbulent_scales

!BOP

! !IROUTINE: cvmix_kpp_compute_shape_function_coeffs
! !INTERFACE:

  subroutine cvmix_kpp_compute_shape_function_coeffs(GAT1, DGAT1, coeffs)

! !DESCRIPTION:
!  Computes the shape function $G(\sigma) = a_0 + a_1\sigma + a_2\sigma^2
!  + a_3\sigma^3$, where
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
    real(cvmix_r8), dimension(4), intent(out) :: coeffs

!EOP
!BOC

    coeffs(1) =  0.0_cvmix_r8
    coeffs(2) =  1.0_cvmix_r8
    coeffs(3) =  3.0_cvmix_r8*GAT1 - DGAT1 - 2.0_cvmix_r8
    coeffs(4) = -2.0_cvmix_r8*GAT1 + DGAT1 + 1.0_cvmix_r8

!EOC

  end subroutine cvmix_kpp_compute_shape_function_coeffs

  function cubic_root_find(coeffs, x0)

    real(cvmix_r8), dimension(4), intent(in) :: coeffs
    real(cvmix_r8),               intent(in) :: x0

    real(cvmix_r8) :: cubic_root_find
    real(cvmix_r8) :: fun_val, root, slope
    integer :: it_cnt

    root = x0
    fun_val = coeffs(1)*(root**3)+coeffs(2)*(root**2)+coeffs(3)*root+coeffs(4)
    do it_cnt = 1, CVMIX_KPP_MAX_NEWTON_ITERS
      if (abs(fun_val).lt.CVMIX_KPP_NEWTON_TOL) &
        exit
      slope = 3.0_cvmix_r8*coeffs(1)*(root**2)+2.0_cvmix_r8*coeffs(2)*root+coeffs(3)
      root = root - fun_val/slope
      fun_val = coeffs(1)*(root**3)+coeffs(2)*(root**2)+coeffs(3)*root+coeffs(4)
    end do
    cubic_root_find = root

  end function cubic_root_find
      
end module cvmix_kpp
