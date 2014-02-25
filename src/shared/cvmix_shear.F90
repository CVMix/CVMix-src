 module cvmix_shear

!BOP
!\newpage
! !MODULE: cvmix_shear
!
! !DESCRIPTION:
!  This module contains routines to initialize the derived types needed for
!  shear mixing, and to set the viscosity and diffusivity coefficients.
!  Presently this scheme has implemented the shear mixing parameterizations
!  from Pacanowski \& Philander (1981) and Large, McWilliams, \& Doney (1994).
!\\
!\\

! !USES:

   use cvmix_kinds_and_types, only : cvmix_r8,                     &
                                     cvmix_zero,                   &
                                     cvmix_one,                    &
                                     cvmix_strlen,                 &
                                     cvmix_data_type
   use cvmix_background, only :      cvmix_bkgnd_params_type,      &
                                     cvmix_bkgnd_lvary_horizontal, &
                                     cvmix_bkgnd_static_diff,      &
                                     cvmix_bkgnd_static_visc
!EOP

   implicit none
   private
   save

!BOP

! !PUBLIC MEMBER FUNCTIONS:

   public :: cvmix_init_shear
   public :: cvmix_coeffs_shear
   public :: cvmix_put_shear
   public :: cvmix_get_shear_real
   public :: cvmix_get_shear_str

   interface cvmix_put_shear
     module procedure cvmix_put_shear_int
     module procedure cvmix_put_shear_real
     module procedure cvmix_put_shear_str
   end interface cvmix_put_shear

! !PUBLIC TYPES:

  ! cvmix_shear_params_type contains the necessary parameters for shear mixing
  ! (currently Pacanowski-Philander or Large et al)
  type, public :: cvmix_shear_params_type
      private
      ! Type of shear mixing to run (PP => Pacanowski-Philander, KPP => LMD94)
      character(len=cvmix_strlen) :: mix_scheme
      ! numerator in viscosity term in PP81
      ! See Eqs. (1) and (2)
      real(cvmix_r8) :: PP_nu_zero  ! units: m^2/s
      ! coefficient of Richardson number in denominator of diff / visc terms
      real(cvmix_r8) :: PP_alpha    ! units: unitless
      ! exponent of denominator in viscosity term
      real(cvmix_r8) :: PP_exp      ! units: unitless
      ! leading coefficient of LMD94 shear mixing formula (max diff / visc)
      ! see Eq. (28b)
      real(cvmix_r8) :: KPP_nu_zero ! units: m^2/s
      ! critical Richardson number value (larger values result in 0 diffusivity
      ! and viscosity)
      real(cvmix_r8) :: KPP_Ri_zero ! units: unitless
      ! Exponent of unitless factor of diff / visc
      real(cvmix_r8) :: KPP_exp     ! units: unitless
  end type cvmix_shear_params_type
!EOP

  type(cvmix_shear_params_type), target :: CVmix_shear_params_saved

 contains

!BOP

! !IROUTINE: cvmix_init_shear
! !INTERFACE:

  subroutine cvmix_init_shear(CVmix_shear_params_user, mix_scheme,            &
                              PP_nu_zero, PP_alpha, PP_exp, KPP_nu_zero,      &
                              KPP_Ri_zero, KPP_exp)

! !DESCRIPTION:
!  Initialization routine for shear (Richardson number-based) mixing. There are
!  currently two supported schemes - set \verb|mix_scheme = 'PP'| to use the
!  Pacanowski-Philander mixing scheme or set \verb|mix_scheme = 'KPP'| to use
!  the interior mixing scheme laid out in Large et al.
!\\
!\\
!  PP requires setting $\nu_0$ (\verb|PP_nu_zero| in this routine), $alpha$ 
!  (\verb|PP_alpha|), and $n$ (\verb|PP_exp|), and returns
!  \begin{eqnarray*}
!  \nu_{PP} & = & \frac{\nu_0}{(1+\alpha \textrm{Ri})^n} + \nu_b \\
!  \kappa_{PP} & = & \frac{\nu}{1+\alpha \textrm{Ri}} + \kappa_b
!  \end{eqnarray*}
!  Note that $\nu_b$ and $\kappa_b$ are set in \verb|cvmix_init_bkgnd()|, which
!  needs to be called separately from this routine.
! \\
! \\
! KPP requires setting $\nu^0$ (\verb|KPP_nu_zero|, $\textrm{Ri}_0 
! ($\verb|KPP_Ri_zero|), and $p_1$ (\verb|KPP_exp|),  and returns
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
    character(len=*), optional, intent(in) :: mix_scheme 
    real(cvmix_r8),   optional, intent(in) :: PP_nu_zero,                     &
                                              PP_alpha,                       &
                                              PP_exp,                         &
                                              KPP_nu_zero,                    &
                                              KPP_Ri_zero,                    &
                                              KPP_exp

! !OUTPUT PARAMETERS:
    type(cvmix_shear_params_type), optional, target, intent(inout) ::         &
                                              CVmix_shear_params_user

!EOP
!BOC

    type(cvmix_shear_params_type), pointer :: CVmix_shear_params_out

    if (present(CVmix_shear_params_user)) then
      CVmix_shear_params_out => CVmix_shear_params_user
    else
      CVmix_shear_params_out => CVmix_shear_params_saved
    end if

    if (present(mix_scheme)) then
      call cvmix_put_shear("mix_scheme", trim(mix_scheme),                    &
                           CVmix_shear_params_user)
    else
      call cvmix_put_shear("mix_scheme", "KPP", CVmix_shear_params_user)
    end if

    select case (trim(CVmix_shear_params_out%mix_scheme))
      case ('PP')
        if (present(PP_nu_zero)) then
          call cvmix_put_shear("PP_nu_zero", PP_nu_zero,                      &
                               CVmix_shear_params_user)
        else
          call cvmix_put_shear("PP_nu_zero", 0.01_cvmix_r8,                   &
                               CVmix_shear_params_user)
        end if

        if (present(PP_alpha)) then
          call cvmix_put_shear("PP_alpha", PP_alpha, CVmix_shear_params_user)
        else
          call cvmix_put_shear("PP_alpha", 5, CVmix_shear_params_user)
        end if

        if (present(PP_exp)) then
          call cvmix_put_shear("PP_exp", PP_exp, CVmix_shear_params_user)
        else
          call cvmix_put_shear("PP_exp", 2, CVmix_shear_params_user)
        end if

      case ('KPP')
        if (present(KPP_nu_zero)) then
          call cvmix_put_shear("KPP_nu_zero", KPP_nu_zero,                    &
                               CVmix_shear_params_user)
        else
          call cvmix_put_shear("KPP_nu_zero", 50e-4_cvmix_r8,                 &
                               CVmix_shear_params_user)
        end if

        if (present(KPP_Ri_zero)) then
          call cvmix_put_shear("KPP_Ri_zero", KPP_Ri_zero,                    &
                               CVmix_shear_params_user)
        else
          call cvmix_put_shear("KPP_Ri_zero", 0.7_cvmix_r8,                   &
                               CVmix_shear_params_user)
        end if

        if (present(KPP_exp)) then
          call cvmix_put_shear("KPP_exp", KPP_exp, CVmix_shear_params_user)
        else
          call cvmix_put_shear("KPP_exp", 3, CVmix_shear_params_user)
        end if

      case DEFAULT
        print*, "ERROR: ", trim(CVmix_shear_params_out%mix_scheme),           &
                " is not a valid choice for shear mixing."
        stop 1

    end select

!EOC

  end subroutine cvmix_init_shear

!BOP

! !IROUTINE: cvmix_coeffs_shear
! !INTERFACE:

  subroutine cvmix_coeffs_shear(CVmix_vars, CVmix_bkgnd_params, colid,        &
                                no_diff, CVmix_shear_params_user)

! !DESCRIPTION:
!  Computes vertical tracer and velocity mixing coefficients for
!  shear-type mixing parameterizations. Note that Richardson number
!  is needed at both T-points and U-points.
!\\
!\\
!
! !USES:
!  only those used by entire module.

! !INPUT PARAMETERS:
    type(cvmix_shear_params_type), target, optional, intent(in) ::            &
                                           CVmix_shear_params_user
    ! PP mixing requires CVmix_bkgnd_params
    type(cvmix_bkgnd_params_type), optional, intent(in) :: CVmix_bkgnd_params
    ! colid is only needed if CVmix_bkgnd_params%lvary_horizontal is true
    integer,                       optional, intent(in) :: colid
    logical,                       optional, intent(in) :: no_diff

! !INPUT/OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars
!EOP
!BOC

    integer                   :: kw ! vertical cell index
    logical                   :: calc_diff
    real(cvmix_r8)            :: nu
    real(cvmix_r8)            :: nu_zero, PP_alpha, KPP_Ri_zero, loc_exp
    real(cvmix_r8)            :: bkgnd_diff, bkgnd_visc
    real(cvmix_r8), pointer, dimension(:) :: RICH
    type(cvmix_shear_params_type), pointer :: CVmix_shear_params

    if (present(CVmix_shear_params_user)) then
      CVmix_shear_params => CVmix_shear_params_user
    else
      CVmix_shear_params => CVmix_shear_params_saved
    end if

    ! Pointer to make the code more legible
    RICH => CVmix_vars%ShearRichardson_iface
    if (.not.present(no_diff)) then
      calc_diff = .true.
    else
      calc_diff = .not.no_diff
    end if

    select case (trim(CVmix_shear_params%mix_scheme))
      case ('PP')
        ! Error checks
        if (.not.present(CVmix_bkgnd_params)) then
          print*, "ERROR: can not run PP mixing without background mixing."
          stop 1
        end if
        if (cvmix_bkgnd_lvary_horizontal(CVmix_bkgnd_params).and.             &
            (.not.present(colid))) then
          print*, "ERROR: background visc and diff vary in horizontal so you",&
                  "must pass column index to cvmix_coeffs_shear"
          stop 1
        end if

        ! Copy parameters to make the code more legible
        nu_zero  = CVmix_shear_params%PP_nu_zero
        PP_alpha = CVmix_shear_params%PP_alpha
        loc_exp  = CVmix_shear_params%PP_exp

        ! Pacanowski-Philander
        do kw=1,CVmix_vars%nlev+1
          bkgnd_diff = cvmix_bkgnd_static_diff(CVmix_bkgnd_params, kw, colid)
          bkgnd_visc = cvmix_bkgnd_static_visc(CVmix_bkgnd_params, kw, colid)
          nu = nu_zero/((cvmix_one+PP_alpha*RICH(kw))**loc_exp)+bkgnd_visc
          CVmix_vars%Mdiff_iface(kw) = nu
          if (calc_diff) &
            CVmix_vars%Tdiff_iface(kw) = nu/(cvmix_one+PP_alpha*RICH(kw)) +   &
                                         bkgnd_diff
        end do

      case ('KPP')
        ! Copy parameters to make the code more legible
        nu_zero     = CVmix_shear_params%KPP_nu_zero
        KPP_Ri_zero = CVmix_shear_params%KPP_Ri_zero
        loc_exp     = CVmix_shear_params%KPP_exp

        ! Large, et al
        do kw=1,CVmix_vars%nlev+1
            if (RICH(kw).lt.0) then
              CVmix_vars%Tdiff_iface(kw) = nu_zero
            else if (RICH(kw).lt.KPP_Ri_zero) then
              CVmix_vars%Tdiff_iface(kw) = nu_zero * (cvmix_one -             &
                   (RICH(kw)/KPP_Ri_zero)**2)**loc_exp
            else ! Ri_g >= Ri_zero
              CVmix_vars%Tdiff_iface(kw) = 0
            end if
        end do
        ! to do: include global params for prandtl number!
        CVmix_vars%Mdiff_iface = CVmix_vars%Tdiff_iface

      case DEFAULT
        ! Note: this error should be caught in cvmix_init_shear
        print*, "ERROR: invalid choice for type of shear mixing."
        stop 1

    end select

!EOC

  end subroutine cvmix_coeffs_shear

!BOP

! !IROUTINE: cvmix_put_shear_int
! !INTERFACE:

  subroutine cvmix_put_shear_int(varname, val, CVmix_shear_params_user)

! !DESCRIPTION:
!  Write an integer value into a cvmix\_shear\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    integer,          intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_shear_params_type), optional, target, intent(inout) ::         &
                                              CVmix_shear_params_user

!EOP
!BOC

    call cvmix_put_shear(varname, real(val,cvmix_r8), CVmix_shear_params_user)

!EOC

  end subroutine cvmix_put_shear_int

!BOP

! !IROUTINE: cvmix_put_shear_real
! !INTERFACE:

  subroutine cvmix_put_shear_real(varname, val, CVmix_shear_params_user)

! !DESCRIPTION:
!  Write a real value into a cvmix\_shear\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    real(cvmix_r8),   intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_shear_params_type), optional, target, intent(inout) ::         &
                                              CVmix_shear_params_user

!EOP
!BOC

    type(cvmix_shear_params_type), pointer :: CVmix_shear_params_out

    if (present(CVmix_shear_params_user)) then
      CVmix_shear_params_out => CVmix_shear_params_user
    else
      CVmix_shear_params_out => CVmix_shear_params_saved
    end if

    select case (trim(varname))
      case ('PP_nu_zero')
        CVmix_shear_params_out%PP_nu_zero = val
      case ('PP_alpha')
        CVmix_shear_params_out%PP_alpha = val
      case ('PP_exp')
        CVmix_shear_params_out%PP_exp = val
      case ('KPP_nu_zero')
        CVmix_shear_params_out%KPP_nu_zero = val
      case ('KPP_Ri_zero')
        CVmix_shear_params_out%KPP_Ri_zero = val
      case ('KPP_exp')
        CVmix_shear_params_out%KPP_exp = val
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end subroutine cvmix_put_shear_real

!BOP

! !IROUTINE: cvmix_put_shear_str
! !INTERFACE:

  subroutine cvmix_put_shear_str(varname, val, CVmix_shear_params_user)

! !DESCRIPTION:
!  Write a string into a cvmix\_shear\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_shear_params_type), optional, target, intent(inout) ::         &
                                              CVmix_shear_params_user

!EOP
!BOC

    type(cvmix_shear_params_type), pointer :: CVmix_shear_params_out

    if (present(CVmix_shear_params_user)) then
      CVmix_shear_params_out => CVmix_shear_params_user
    else
      CVmix_shear_params_out => CVmix_shear_params_saved
    end if

    select case (trim(varname))
      case ('mix_scheme')
        CVmix_shear_params_out%mix_scheme = val
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end subroutine cvmix_put_shear_str

!BOP

! !IROUTINE: cvmix_get_shear_real
! !INTERFACE:

  function cvmix_get_shear_real(varname, CVmix_shear_params_user)

! !DESCRIPTION:
!  Read the real value of a cvmix\_shear\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*),                                intent(in) :: varname
    type(cvmix_shear_params_type), optional, target, intent(in) ::            &
                                           CVmix_shear_params_user

! !OUTPUT PARAMETERS:
    real(cvmix_r8) :: cvmix_get_shear_real

!EOP
!BOC

    type(cvmix_shear_params_type), pointer :: CVmix_shear_params_in

    if (present(CVmix_shear_params_user)) then
      CVmix_shear_params_in => CVmix_shear_params_user
    else
      CVmix_shear_params_in => CVmix_shear_params_saved
    end if

    cvmix_get_shear_real = cvmix_zero
    select case (trim(varname))
      case ('PP_nu_zero')
        cvmix_get_shear_real =CVmix_shear_params_in%PP_nu_zero
      case ('PP_alpha')
        cvmix_get_shear_real =CVmix_shear_params_in%PP_alpha
      case ('PP_exp')
        cvmix_get_shear_real =CVmix_shear_params_in%PP_exp
      case ('KPP_nu_zero')
        cvmix_get_shear_real =CVmix_shear_params_in%KPP_nu_zero
      case ('KPP_Ri_zero')
        cvmix_get_shear_real =CVmix_shear_params_in%KPP_Ri_zero
      case ('KPP_exp')
        cvmix_get_shear_real =CVmix_shear_params_in%KPP_exp
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end function cvmix_get_shear_real

!BOP

! !IROUTINE: cvmix_get_shear_str
! !INTERFACE:

  function cvmix_get_shear_str(varname, CVmix_shear_params_user)

! !DESCRIPTION:
!  Read the string contents of a cvmix\_shear\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*),                                intent(in) :: varname
    type(cvmix_shear_params_type), optional, target, intent(in) ::            &
                                           CVmix_shear_params_user

! !OUTPUT PARAMETERS:
    character(len=cvmix_strlen) :: cvmix_get_shear_str

!EOP
!BOC

    type(cvmix_shear_params_type), pointer :: CVmix_shear_params_in

    if (present(CVmix_shear_params_user)) then
      CVmix_shear_params_in => CVmix_shear_params_user
    else
      CVmix_shear_params_in => CVmix_shear_params_saved
    end if

    select case (trim(varname))
      case ('mix_scheme')
        cvmix_get_shear_str = trim(CVmix_shear_params_in%mix_scheme)
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end function cvmix_get_shear_str


end module cvmix_shear
