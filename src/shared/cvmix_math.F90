!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module cvmix_math

!BOP
!\newpage
! !MODULE: cvmix_math
!
! !DESCRIPTION:
!  This module contains routines to compute polynomial interpolations (linear, 
!  quadratic, or cubic spline), evaluate  third-order polynomials and their
!  derivatives at specific values, and compute roots of these polynomials.
!\\
!\\
!
! !REVISION HISTORY:
!  $Id$
!  $URL$

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8

!EOP

  implicit none
  private
  save

!BOP

! !DEFINED PARAMETERS:
  integer, parameter, public :: CVMIX_MATH_INTERP_LINEAR      = 1
  integer, parameter, public :: CVMIX_MATH_INTERP_QUAD        = 2
  integer, parameter, public :: CVMIX_MATH_INTERP_CUBE_SPLINE = 3

  real(cvmix_r8), parameter :: CVMIX_MATH_NEWTON_TOL       = 1.0e-12_cvmix_r8
  integer,        parameter :: CVMIX_MATH_MAX_NEWTON_ITERS = 100

! !PUBLIC MEMBER FUNCTIONS:

  public :: cvmix_math_poly_interp
  public :: cvmix_math_cubic_root_find
  public :: cvmix_math_evaluate_cubic

!EOP

  contains

!BOP

! !IROUTINE: cvmix_math_poly_interp
! !INTERFACE:

  subroutine cvmix_math_poly_interp(coeffs, interp_type, x, y, x0, y0)

! !INPUT PARAMETERS:
    integer,                      intent(in)    :: interp_type
    real(cvmix_r8), dimension(2), intent(in)    :: x, y
    real(cvmix_r8), optional,     intent(in)    :: x0, y0
! !OUTPUT PARAMETERS:
    real(cvmix_r8), dimension(4), intent(inout) :: coeffs

!EOP
!BOC

    ! Local variables
    real(cvmix_r8) :: det
    integer        :: k, k2
    real(kind=cvmix_r8), dimension(:,:), allocatable :: Minv
    real(kind=cvmix_r8), dimension(:),   allocatable :: rhs

    ! All interpolation assumes form of
    ! y = ax^3 + bx^2 + cx + d
    ! linear => a = b = 0
    ! quad   => a = 0
    coeffs(1:4) = 0.0_cvmix_r8
    select case (interp_type)
      case (CVMIX_MATH_INTERP_LINEAR)
        ! Match y(1) and y(2)
!        print*, "Linear interpolation"
        coeffs(3) = (y(2)-y(1))/(x(2)-x(1))
        coeffs(4) = y(1)-coeffs(3)*x(1)
      case (CVMIX_MATH_INTERP_QUAD)
        ! Match y(1), y(2), and y'(1) [requires x(0)]
!        print*, "Quadratic interpolation"
        ! [ x2^2 x2 1 ][ b ]   [    y2 ]
        ! [ x1^2 x1 1 ][ c ] = [    y1 ]
        ! [  2x1  1 0 ][ d ]   [ slope ]
        !      ^^^
        !       M
        det = -((x(2)-x(1))**2)
        allocate(Minv(3,3))
        allocate(rhs(3))
        rhs(1) = y(2)
        rhs(2) = y(1)
        if (present(x0).and.present(y0)) then
          rhs(3) = (y(1)-y0)/(x(1)-x0)
        else
          rhs(3) = 0.0_cvmix_r8
        end if

        Minv(1,1) = -real(1, cvmix_r8)/det
        Minv(1,2) = real(1, cvmix_r8)/det
        Minv(1,3) = -real(1, cvmix_r8)/(x(2)-x(1))
        Minv(2,1) = real(2, cvmix_r8)*x(1)/det
        Minv(2,2) = -real(2, cvmix_r8)*x(1)/det
        Minv(2,3) = (x(2)+x(1))/(x(2)-x(1))
        Minv(3,1) = -(x(1)**2)/det
        Minv(3,2) = x(2)*(real(2, cvmix_r8)*x(1)-x(2))/det
        Minv(3,3) = -x(2)*x(1)/(x(2)-x(1))

        do k=1,3
          coeffs(2) = coeffs(2)+Minv(1,k)*rhs(k)
          coeffs(3) = coeffs(3)+Minv(2,k)*rhs(k)
          coeffs(4) = coeffs(4)+Minv(3,k)*rhs(k)
        end do
        deallocate(rhs)
        deallocate(Minv)
      case (CVMIX_MATH_INTERP_CUBE_SPLINE)
        ! Match y(1), y(2), y'(1), and y'(2)
!        print*, "Cubic spline interpolation"
        ! [ x2^3 x2^2 x2 1 ][ a ]   [     y2 ]
        ! [ x1^3 x1^2 x1 1 ][ b ] = [     y1 ]
        ! [  3x1  2x1  1 0 ][ c ]   [ slope1 ]
        ! [  3x2  2x2  1 0 ][ d ]   [ slope2 ]
        !      ^^^
        !       M
        det = -((x(2)-x(1))**3)
        allocate(Minv(4,4))
        allocate(rhs(4))
        rhs(1) = y(2)
        rhs(2) = y(1)
        if (present(x0).and.present(y0)) then
          rhs(3) = (y(1)-y0)/(x(1)-x0)
        else
          rhs(3) = 0.0_cvmix_r8
        end if
        rhs(4) = (y(2)-y(1))/(x(2)-x(1))

        Minv(1,1) = real(2, cvmix_r8)/det
        Minv(1,2) = -real(2, cvmix_r8)/det
        Minv(1,3) = (x(1)-x(2))/det
        Minv(1,4) = (x(1)-x(2))/det
        Minv(2,1) = -real(3, cvmix_r8)*(x(2)+x(1))/det
        Minv(2,2) = real(3, cvmix_r8)*(x(2)+x(1))/det
        Minv(2,3) = (x(2)-x(1))*(real(2, cvmix_r8)*x(2)+x(1))/det
        Minv(2,4) = (x(2)-x(1))*(real(2, cvmix_r8)*x(1)+x(2))/det
        Minv(3,1) = real(6, cvmix_r8)*x(2)*x(1)/det
        Minv(3,2) = -real(6, cvmix_r8)*x(2)*x(1)/det
        Minv(3,3) = -x(2)*(x(2)-x(1))*(real(2, cvmix_r8)*x(1)+x(2))/det
        Minv(3,4) = -x(1)*(x(2)-x(1))*(real(2, cvmix_r8)*x(2)+x(1))/det
        Minv(4,1) = -(x(1)**2)*(real(3, cvmix_r8)*x(2)-x(1))/det
        Minv(4,2) = -(x(2)**2)*(-real(3, cvmix_r8)*x(1)+x(2))/det
        Minv(4,3) = x(1)*(x(2)**2)*(x(2)-x(1))/det
        Minv(4,4) = x(2)*(x(1)**2)*(x(2)-x(1))/det

        do k=1,4
          do k2=1,4
            coeffs(k2) = coeffs(k2)+Minv(k2,k)*rhs(k)
          end do
        end do
        deallocate(rhs)
        deallocate(Minv)
    end select

!EOC

  end subroutine cvmix_math_poly_interp

  function cvmix_math_cubic_root_find(coeffs, x0)

    real(cvmix_r8), dimension(4), intent(in) :: coeffs
    real(cvmix_r8),               intent(in) :: x0

    real(cvmix_r8) :: cvmix_math_cubic_root_find
    real(cvmix_r8) :: fun_val, root, slope
    integer :: it_cnt

    root = x0
    fun_val = coeffs(1)*(root**3)+coeffs(2)*(root**2)+coeffs(3)*root+coeffs(4)
    do it_cnt = 1, CVMIX_MATH_MAX_NEWTON_ITERS
      if (abs(fun_val).lt.CVMIX_MATH_NEWTON_TOL) &
        exit
      slope = 3.0_cvmix_r8*coeffs(1)*(root**2)+2.0_cvmix_r8*coeffs(2)*root+coeffs(3)
      root = root - fun_val/slope
      fun_val = coeffs(1)*(root**3)+coeffs(2)*(root**2)+coeffs(3)*root+coeffs(4)
    end do
    cvmix_math_cubic_root_find = root

  end function cvmix_math_cubic_root_find
      
!BOP

! !IROUTINE: cvmix_math_evaluate_cubic
! !INTERFACE:

  function cvmix_math_evaluate_cubic(coeffs, x_in, fprime)

! !DESCRIPTION:
!  Computes $f(x) = a_0 + a_1x + a_2x^2 + a_3x^3$ at $x = $\verb|x_in|, where
!  \verb|coeffs|$ = (/a_0, a_1, a_2, a_3/)$. If requested, can also return
!  $f'(x)$
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    real(cvmix_r8), dimension(4), intent(in) :: coeffs
    real(cvmix_r8),               intent(in) :: x_in

! !OUTPUT PARAMETERS:
    real(cvmix_r8) :: cvmix_math_evaluate_cubic
    real(cvmix_r8), optional, intent(out) :: fprime

!EOP
!BOC

    ! Local Variables
    integer :: i

    cvmix_math_evaluate_cubic = 0.0_cvmix_r8
      if (present(fprime)) &
        fprime = 0.0_cvmix_r8
    do i=1,4
      cvmix_math_evaluate_cubic = cvmix_math_evaluate_cubic +                 &
                                  coeffs(i)*(x_in**(i-1))
      if (present(fprime).and.(i.gt.1)) &
        fprime = fprime + real(i-1,cvmix_r8)*(x_in**(i-2))
    end do

  end function cvmix_math_evaluate_cubic

end module cvmix_math
