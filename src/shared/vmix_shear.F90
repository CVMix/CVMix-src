!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module vmix_shear

!BOP
!\newpage
! !MODULE: vmix_shear
!
! !DESCRIPTION:
!  This module contains routines to initialize the derived types needed for
!  shear mixing (currently just the Pacanowski-Philander scheme) and to set
!  the viscosity and diffusivity coefficients accordingly.
!\\
!\\
!
! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !USES:

   use vmix_kinds_and_types, only : vmix_r8,                  &
                                    vmix_strlen,              &
                                    vmix_data_type,           &
                                    vmix_global_params_type,  &
                                    vmix_bkgnd_params_type,   &
                                    vmix_shear_params_type
!EOP

   implicit none
   private
   save

!BOP

! !PUBLIC MEMBER FUNCTIONS:

   public :: vmix_init_shear
   public :: vmix_coeffs_shear
!EOP

 contains

!BOP

! !IROUTINE: vmix_init_shear
! !INTERFACE:

  subroutine vmix_init_shear(Vmix_shear_params, mix_scheme,     &
                             alpha, n, nu_zero, Ri_zero, p_one)

! !DESCRIPTION:
!  Initialization routine for shear (Richardson number-based) mixing. There are
!  currently two supported schemes - set \verb|mix_scheme = 'PP'| to use the
!  Paconowski-Philander mixing scheme or set \verb|mix_scheme = 'KPP'| to use
!  the interior mixing scheme laid out in Large et al.
!\\
!\\
!  PP requires setting \verb|nu_zero| ($\nu_0$), \verb|alpha| ($\alpha$), and
!  \verb|n| ($n$), and returns
!  \begin{eqnarray*}
!  \nu_{PP} & = & \frac{\nu_0}{(1+\alpha \textrm{Ri})^n} + \nu_b \\
!  \kappa_{PP} & = & \frac{\nu}{1+\alpha \textrm{Ri}} + \kappa_b
!  \end{eqnarray*}
!  Note that $\nu_b$ and $\kappa_b$ are set in \verb|vmix_init_bkgnd()|, which
!  needs to be called separately from this routine.
! \\
! \\
! KPP requires setting \verb|nu_zero| ($\nu^0$), \verb|p_one| ($p_1$), and
! \verb|Ri_zero| ($\textrm{Ri}_0$), and returns
! $$
! \nu_{KPP} = \left\{
! \begin{array}{r l}
! \nu^0 & \textrm{Ri} < 0\\
! \nu^0 \left[1 - \frac{\textrm{Ri}}{\textrm{Ri}_0}^2\right]^{p_1}
!       & 0 < \textrm{Ri}
!           < \textrm{Ri}_0 \\
! 0     & \textrm{Ri}_0 < \textrm{Ri}
! \end{array} \right.
! $$
!
! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*),        intent(in) :: mix_scheme 
    real(vmix_r8), optional, intent(in) :: alpha, n, nu_zero, ri_zero, p_one

! !OUTPUT PARAMETERS:
    type(vmix_shear_params_type), intent(inout) :: Vmix_shear_params
!EOP
!BOC

    select case (trim(mix_scheme))
      case ('PP')
        if (.not.(present(alpha).and.present(nu_zero).and.present(n))) then
          print*, "ERROR: you must specify alpha, nu_zero, and n to use" ,&
                  "Paconowski-Philander mixing!"
          stop
        end if
        Vmix_shear_params%mix_scheme = "PP"
        Vmix_shear_params%alpha      = alpha
        Vmix_shear_params%n          = n
        Vmix_shear_params%nu_zero    = nu_zero

      case ('KPP')
        Vmix_shear_params%mix_scheme = "PP"
        Vmix_shear_params%nu_zero    = nu_zero
        Vmix_shear_params%Ri_zero    = Ri_zero
        Vmix_shear_params%p_one      = p_one

      case DEFAULT
        print*, "ERROR: ", trim(mix_scheme), " is not a valid choice for ", &
                "shear mixing."
        stop

    end select

!EOC

  end subroutine vmix_init_shear

!***********************************************************************
!BOP
! !IROUTINE: vmix_coeffs_shear
! !INTERFACE:

  subroutine vmix_coeffs_shear(Vmix_vars, Vmix_shear_params, &
                               Vmix_bkgnd_params, colid)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for shear-type mixing
!  parameterizatiions. Note that Richardson number is needed at
!  both T-points and U-points.
!\\
!\\
!
! !USES:
!  only those used by entire module.

! !INPUT PARAMETERS:
    type(vmix_shear_params_type), intent(in) :: Vmix_shear_params
    type(vmix_bkgnd_params_type), intent(in) :: Vmix_bkgnd_params
    ! colid is only needed if Vmix_bkgnd_params%lvary_horizontal is true
    integer, optional,            intent(in) :: colid

! !INPUT/OUTPUT PARAMETERS:
    type(vmix_data_type), intent(inout) :: Vmix_vars
!EOP
!BOC

    integer                  :: kw ! vertical cell index
    real(vmix_r8), parameter :: one = 1.0_vmix_r8
    real(vmix_r8)            :: nu
    real(vmix_r8)            :: alpha, n, nu_zero, Ri_zero, p_one
    real(vmix_r8)            :: bkgnd_diff, bkgnd_visc
    real(vmix_r8), pointer, dimension(:) :: RICHT, RICHU

    ! Copy vars / create pointers to make the code more legible
    alpha   = Vmix_shear_params%alpha
    n       = Vmix_shear_params%n
    nu_zero = Vmix_shear_params%nu_zero
    Ri_zero = Vmix_shear_params%Ri_zero
    p_one   = Vmix_shear_params%p_one
    RICHT => Vmix_vars%Ri_t_iface
    RICHU => Vmix_vars%Ri_u_iface

    ! Error check
    if (Vmix_bkgnd_params%lvary_horizontal.and.(.not.present(colid))) then
      print*, "ERROR: background visc and diff vary in horizontal so you", &
              "must pass column index to vmix_coeffs_shear"
      stop
    end if

    select case (trim(Vmix_shear_params%mix_scheme))
      case ('PP')
        ! Paconowski-Philander
        do kw=1,Vmix_vars%nlev+1
          if (Vmix_bkgnd_params%lvary_horizontal) then
            if (Vmix_bkgnd_params%lvary_vertical) then
              bkgnd_diff = Vmix_bkgnd_params%static_diff(colid, kw)
              bkgnd_visc = Vmix_bkgnd_params%static_visc(colid, kw)
            else
              bkgnd_diff = Vmix_bkgnd_params%static_diff(colid, 1)
              bkgnd_visc = Vmix_bkgnd_params%static_visc(colid, 1)
            end if
          else
            if (Vmix_bkgnd_params%lvary_vertical) then
              bkgnd_diff = Vmix_bkgnd_params%static_diff(1, kw)
              bkgnd_visc = Vmix_bkgnd_params%static_visc(1, kw)
            else
              bkgnd_diff = Vmix_bkgnd_params%static_diff(1, 1)
              bkgnd_visc = Vmix_bkgnd_params%static_visc(1, 1)
            end if
          end if
          nu = nu_zero/((one+alpha*RICHT(kw))**n)+bkgnd_visc
          Vmix_vars%diff_iface(kw,1) = nu/((one+alpha*RICHT(kw))**n) +    &
                                       bkgnd_diff
          Vmix_vars%visc_iface(kw) = nu_zero/((one+alpha*RICHU(kw))**n) + &
                                     bkgnd_visc
        end do

      case ('KPP')
        ! Large, et al
        do kw=1,Vmix_vars%nlev+1
            ! On T-levs only (need to figure out how to split T and U)
            if (RICHT(kw).lt.0) then
              Vmix_vars%diff_iface(kw,1) = nu_zero
            else if (RICHT(kw).lt.Ri_zero) then
              Vmix_vars%diff_iface(kw,1) = nu_zero * (one -               &
                   (RICHT(kw)/Ri_zero)**2)**p_one 
            else ! Ri_g >= Ri_zero
              Vmix_vars%diff_iface(kw,1) = 0
            end if
        end do

      case DEFAULT
        ! Note: this error should be caught in vmix_init_shear
        print*, "ERROR: invalid choice for type of shear mixing."
        stop

    end select

  end subroutine vmix_coeffs_shear

end module vmix_shear
