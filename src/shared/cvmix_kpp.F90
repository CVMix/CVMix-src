!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module cvmix_kpp

!BOP
!\newpage
! !MODULE: cvmix_kpp
!
! !AUTHOR:
!  Michael N. Levy, NCAR (mlevy@ucar.edu)
!
! !DESCRIPTION:
!  This module contains routines to initialize the derived types needed for
!  KPP mixing and to set the viscosity and diffusivity coefficients
!  accordingly.
!\\
!\\
!  References:\\
!  * WG Large, JC McWilliams, and SC Doney.
!  Oceanic Vertical Mixing: A Review and a Model with a Nonlocal Boundary Layer
!  Parameterization.
!  Review of Geophysics, 1994.
!\\
!\\

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                                 &
                                    cvmix_strlen,                             &
                                    cvmix_zero,                               &
                                    cvmix_one,                                &
                                    cvmix_PI,                                 &
                                    cvmix_data_type,                          &
                                    cvmix_global_params_type,                 &
                                    CVMIX_OVERWRITE_OLD_VAL,                  &
                                    CVMIX_SUM_OLD_AND_NEW_VALS,               &
                                    CVMIX_MAX_OLD_AND_NEW_VALS
  use cvmix_math, only :            CVMIX_MATH_INTERP_LINEAR,                 &
                                    CVMIX_MATH_INTERP_QUAD,                   &
                                    CVMIX_MATH_INTERP_CUBE_SPLINE,            &
                                    cvmix_math_poly_interp,                   &
                                    cvmix_math_cubic_root_find,               &
                                    cvmix_math_evaluate_cubic
  use cvmix_put_get,         only : cvmix_put
  use cvmix_utils,           only : cvmix_update_wrap

!EOP

  implicit none
  private
  save

!BOP

! !DEFINED PARAMETERS:
  integer, parameter :: CVMIX_KPP_INTERP_LMD94       = -1
  integer, parameter :: CVMIX_KPP_MATCH_BOTH         = 1
  integer, parameter :: CVMIX_KPP_MATCH_GRADIENT     = 2
  integer, parameter :: CVMIX_KPP_SIMPLE_SHAPES      = 3
  integer, parameter :: CVMIX_KPP_PARABOLIC_NONLOCAL = 4
  integer, parameter :: NO_LANGMUIR_MIXING           = -1
  integer, parameter :: LANGMUIR_MIXING_LWF16        = 1
  integer, parameter :: LANGMUIR_MIXING_RWHGK16      = 2
  integer, parameter :: NO_LANGMUIR_ENTRAINMENT      = -1
  integer, parameter :: LANGMUIR_ENTRAINMENT_LWF16   = 1
  integer, parameter :: LANGMUIR_ENTRAINMENT_LF17    = 2
  integer, parameter :: LANGMUIR_ENTRAINMENT_RWHGK16 = 3

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
  public :: cvmix_kpp_composite_Gshape                     !STOKES_MOST
  public :: cvmix_kpp_compute_StokesXi                    !STOKES_MOST
  ! These are public for testing, may end up private later
  public :: cvmix_kpp_compute_shape_function_coeffs
  public :: cvmix_kpp_compute_kOBL_depth
  public :: cvmix_kpp_compute_enhanced_diff
  public :: cvmix_kpp_compute_nu_at_OBL_depth_LMD94
  public :: cvmix_kpp_EFactor_model
  public :: cvmix_kpp_ustokes_SL_model


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
    module procedure cvmix_kpp_compute_turbulent_scales_1d_sigma
    module procedure cvmix_kpp_compute_turbulent_scales_1d_OBL
  end interface cvmix_kpp_compute_turbulent_scales

! !PUBLIC TYPES:

  ! cvmix_kpp_params_type contains the necessary parameters for KPP mixing
  type, public :: cvmix_kpp_params_type
    private
      real(cvmix_r8) :: Ri_crit        ! Critical Richardson number
                                       ! (OBL_depth = where bulk Ri = Ri_crit)

      real(cvmix_r8) :: minOBLdepth    ! Minimum allowable OBL depth
                                       ! (Default is 0 m => no minimum)
      real(cvmix_r8) :: maxOBLdepth    ! Maximum allowable OBL depth
                                       ! (Default is 0 m => no maximum)
      real(cvmix_r8) :: minVtsqr       ! Minimum allowable unresolved shear
                                       ! (Default is 1e-10 m^2/s^2)

      real(cvmix_r8) :: vonkarman      ! von Karman constant

      real(cvmix_r8) :: Cstar          ! coefficient for nonlinear transport
      real(cvmix_r8) :: nonlocal_coeff ! Cs from Eq (20) in LMD94
                                       ! Default value comes from paper, but
                                       ! some users may set it = 1.

      ! For velocity scale function, _m => momentum and _s => scalar (tracer)
      real(cvmix_r8) :: zeta_m         ! parameter for computing vel scale func
      real(cvmix_r8) :: zeta_s         ! parameter for computing vel scale func
      real(cvmix_r8) :: a_m            ! parameter for computing vel scale func
      real(cvmix_r8) :: c_m            ! parameter for computing vel scale func
      real(cvmix_r8) :: a_s            ! parameter for computing vel scale func
      real(cvmix_r8) :: c_s            ! parameter for computing vel scale func

      real(cvmix_r8) :: surf_layer_ext ! nondimensional extent of surface layer
                                       ! (expressed in sigma-coordinates)

      integer        :: interp_type    ! interpolation type used to interpolate
                                       ! bulk Richardson number
      integer        :: interp_type2   ! interpolation type used to interpolate
                                       ! diff and visc at OBL_depth

      ! Cv is a parameter used to compute the unresolved shear. By default, the
      ! formula from Eq. (A3) of Danabasoglu et al. is used, but a single
      ! scalar value can be set instead.
      real(cvmix_r8) :: Cv

      ! MatchTechnique is set by a string of the same name as an argument in
      ! cvmix_init_kpp. It determines how matching between the boundary layer
      ! and ocean interior is handled at the interface. Note that this also
      ! controls whether the shape function used to compute the coefficient in
      ! front of the nonlocal term is the same as that used to compute the
      ! gradient term.
      ! Options (for cvmix_init_kpp) are
      ! (i) SimpleShapes => Shape functions for both the gradient and nonlocal
      !                     terms vanish at interface
      ! (ii) MatchGradient => Shape function for nonlocal term vanishes at
      !                       interface, but gradient term matches interior
      !                       values.
      ! (iii) MatchBoth => Shape functions for both the gradient and nonlocal
      !                    term match interior values at interface
      ! (iv) ParabolicNonLocal => Shape function for the nonlocal term is
      !                         (1-sigma)^2, gradient term is sigma*(1-sigma)^2
      integer :: MatchTechnique

      ! Flag for what to do with old values of CVmix_vars%[MTS]diff
      integer :: handle_old_vals

      ! Logic flags to dictate if / how various terms are computed
      logical        :: lStokesMOST    ! True => use Stokes Similarty package
      logical        :: lscalar_Cv     ! True => use the scalar Cv value
      logical        :: lEkman         ! True => compute Ekman depth limit
      logical        :: lMonOb         ! True => compute Monin-Obukhov limit
      logical        :: lnoDGat1       ! True => G'(1) = 0 (shape function)
                                       ! False => compute G'(1) as in LMD94
      logical        :: lenhanced_diff ! True => enhance diffusivity at OBL
      integer        :: Langmuir_Mixing_Opt
                                       ! Option of Langmuir enhanced mixing
                                       ! - apply an enhancement factor to the
                                       ! turbulent velocity scale
      integer        :: Langmuir_Entrainment_Opt
                                       ! Option of Langmuir turbulence enhanced
                                       ! entrainment - modify the unresolved shear
      logical        :: l_LMD_ws       ! flag to use original Large et al. (1994)
                                       ! equations for computing turbulent scales
                                       ! rather than the updated methodology in
                                       ! Danabasoglu et al. (2006). The latter
                                       ! limits sigma to be < surf_layer_extent
                                       ! when computing turbulent scales while
                                       ! the former only imposes this restriction
                                       ! in unstable regimes.
      real(cvmix_r8) :: c_LT, c_ST, c_CT  ! Empirical constants in the scaling of the
                                          ! entrainment buoyancy flux
                                          ! (20) in Li and Fox-Kemper, 2017, JPO
      real(cvmix_r8) :: p_LT              ! Power of Langmuir number in the above
                                          ! scaling
      !BGR
      real(cvmix_r8) :: RWHGK_ENTR_COEF,& ! Coefficient and exponent from
                        RWHGK_ENTR_EXP    ! RWHGK16 Langmuir parameterization

      real(cvmix_r8) :: CVt2              ! Tunable parameter for convection entrainment
                                          ! (Only used with StokesMOST)

  end type cvmix_kpp_params_type

!EOP

type(cvmix_kpp_params_type), target :: CVmix_kpp_params_saved

contains

!BOP

! !IROUTINE: cvmix_init_kpp
! !INTERFACE:
  subroutine cvmix_init_kpp(ri_crit, minOBLdepth, maxOBLdepth, minVtsqr,      &
                            vonkarman, Cstar, zeta_m, zeta_s, surf_layer_ext, &
                            Cv, interp_type, interp_type2, MatchTechnique,    &
                            old_vals, lEkman, lStokesMOST, lMonOb, lnoDGat1,  &
                            lenhanced_diff, lnonzero_surf_nonlocal,           &
                            Langmuir_mixing_str, Langmuir_entrainment_str,    &
                            l_LMD_ws, CVmix_kpp_params_user)

! !DESCRIPTION:
!  Initialization routine for KPP mixing.
!\\
!\\

! !INPUT PARAMETERS:
    real(cvmix_r8),   optional, intent(in) :: ri_crit,                        &
                                              minOBLdepth,                    &
                                              maxOBLdepth,                    &
                                              minVtsqr,                       &
                                              vonkarman,                      &
                                              Cstar,                          &
                                              zeta_m,                         &
                                              zeta_s,                         &
                                              surf_layer_ext,                 &
                                              Cv
    character(len=*), optional, intent(in) :: interp_type,                    &
                                              interp_type2,                   &
                                              MatchTechnique,                 &
                                              old_vals,                       &
                                              Langmuir_mixing_str,            &
                                              Langmuir_entrainment_str
    logical,          optional, intent(in) :: lEkman,                         &
                                              lStokesMOST,                    &
                                              lMonOb,                         &
                                              lnoDGat1,                       &
                                              lenhanced_diff,                 &
                                              lnonzero_surf_nonlocal,         &
                                              l_LMD_ws

! !OUTPUT PARAMETERS:
    type(cvmix_kpp_params_type), intent(inout), target, optional ::           &
                                              CVmix_kpp_params_user

!EOP
!BOC

    real(cvmix_r8) :: zm, zs, a_m, a_s, c_m, c_s
    real(cvmix_r8) :: Cstar_loc, vonkar_loc, surf_layer_ext_loc
    real(cvmix_r8) :: nonlocal_coeff

    if (present(ri_crit)) then
      if (ri_crit.lt.cvmix_zero) then
        print*, "ERROR: ri_crit can not be negative."
        stop 1
      end if
      call cvmix_put_kpp('Ri_crit', ri_crit, CVmix_kpp_params_user)
    else
      call cvmix_put_kpp('Ri_crit', 0.3_cvmix_r8, CVmix_kpp_params_user)
    end if

    if (present(minOBLdepth)) then
      if (minOBLdepth.lt.cvmix_zero) then
        print*, "ERROR: minOBLdepth can not be negative."
        stop 1
      end if
      call cvmix_put_kpp('minOBLdepth', minOBLdepth, CVmix_kpp_params_user)
    else
      call cvmix_put_kpp('minOBLdepth', 0, CVmix_kpp_params_user)
    end if

    if (present(maxOBLdepth)) then
      if (maxOBLdepth.lt.cvmix_zero) then
        print*, "ERROR: maxOBLdepth can not be negative."
        stop 1
      end if
      call cvmix_put_kpp('maxOBLdepth', maxOBLdepth, CVmix_kpp_params_user)
    else
      call cvmix_put_kpp('maxOBLdepth', 0, CVmix_kpp_params_user)
    end if

    if (present(minVtsqr)) then
      if (minVtsqr.lt.cvmix_zero) then
        print*, "ERROR: minVtsqr can not be negative."
        stop 1
      end if
      call cvmix_put_kpp('minVtsqr', minVtsqr, CVmix_kpp_params_user)
    else
      call cvmix_put_kpp('minVtsqr', 1e-10_cvmix_r8, CVmix_kpp_params_user)
    end if

    if (present(vonkarman)) then
      if (vonkarman.lt.cvmix_zero) then
        print*, "ERROR: vonkarman can not be negative."
        stop 1
      end if
      vonkar_loc = vonkarman
    else
      vonkar_loc = 0.4_cvmix_r8
    end if
    call cvmix_put_kpp('vonkarman', vonkar_loc, CVmix_kpp_params_user)

    if (present(Cstar)) then
      Cstar_loc = Cstar
    else
      Cstar_loc = real(10,cvmix_r8)
    end if
    call cvmix_put_kpp('Cstar', Cstar_loc, CVmix_kpp_params_user)

    if (present(zeta_m)) then
      if (zeta_m.ge.cvmix_zero) then
        print*, "ERROR: zeta_m must be negative."
        stop 1
      end if
      zm = zeta_m
    else
      ! default value for zeta_m is -1/5
      zm = -0.2_cvmix_r8
    end if
    call cvmix_put_kpp('zeta_m', zm, CVmix_kpp_params_user)

    if (present(zeta_s)) then
      if (zeta_s.ge.cvmix_zero) then
        print*, "ERROR: zeta_s must be negative."
        stop 1
      end if
      zs = zeta_s
    else
      ! Default value for zeta_s is -1
      zs = -cvmix_one
    end if
    call cvmix_put_kpp('zeta_s', zs, CVmix_kpp_params_user)

    ! a_m, a_s, c_m, and c_s are computed from zeta_m and zeta_s
    ! a_m, c_m, and c_s are all non-negative. a_s may be negative depending
    ! on the value of zeta_s
    a_m = ((cvmix_one - real(16,cvmix_r8)*zm)**(-0.25_cvmix_r8))*             &
          (cvmix_one - real(4,cvmix_r8)*zm)
    c_m = ((cvmix_one - real(16,cvmix_r8)*zm)**(-0.25_cvmix_r8))*             &
          real(12,cvmix_r8)
    call cvmix_put_kpp('a_m', a_m, CVmix_kpp_params_user)
    call cvmix_put_kpp('c_m', c_m, CVmix_kpp_params_user)

    a_s = sqrt(cvmix_one - real(16,cvmix_r8)*zs)*                             &
          (cvmix_one + real(8,cvmix_r8)*zs)
    c_s = real(24,cvmix_r8)*sqrt(cvmix_one - real(16,cvmix_r8)*zs)
    call cvmix_put_kpp('a_s', a_s, CVmix_kpp_params_user)
    call cvmix_put_kpp('c_s', c_s, CVmix_kpp_params_user)

    if (present(surf_layer_ext)) then
      if ((surf_layer_ext.lt.cvmix_zero).or.(surf_layer_ext.gt.cvmix_one))    &
      then
        print*, "surf_layer_ext must be between 0 and 1, inclusive."
        stop 1
      end if
      surf_layer_ext_loc = surf_layer_ext
    else
      surf_layer_ext_loc = 0.1_cvmix_r8
    end if
    call cvmix_put_kpp('surf_layer_ext', surf_layer_ext_loc,                  &
                       CVmix_kpp_params_user)

    if (present(Cv)) then
      ! Use scalar Cv parameter
      call cvmix_put_kpp('Cv', CV, CVmix_kpp_params_user)
      call cvmix_put_kpp('lscalar_Cv', .true., CVmix_kpp_params_user)
    else
      ! Use Eq. (A3) from Danabasoglu et al.
      call cvmix_put_kpp('lscalar_Cv', .false., CVmix_kpp_params_user)
    end if

    if (present(interp_type)) then
      select case (trim(interp_type))
        case ('line', 'linear')
          call cvmix_put_kpp('interp_type', CVMIX_MATH_INTERP_LINEAR,         &
                             CVmix_kpp_params_user)
        case ('quad', 'quadratic')
          call cvmix_put_kpp('interp_type', CVMIX_MATH_INTERP_QUAD,           &
                             CVmix_kpp_params_user)
        case ('cube', 'cubic', 'cubic_spline', 'cubic spline')
          call cvmix_put_kpp('interp_type', CVMIX_MATH_INTERP_CUBE_SPLINE,    &
                             CVmix_kpp_params_user)
        case DEFAULT
          print*, "ERROR: ", trim(interp_type), " is not a valid type of ",   &
                  "interpolation!"
          stop 1
      end select
    else
      call cvmix_put_kpp('interp_type', CVMIX_MATH_INTERP_QUAD,               &
                         CVmix_kpp_params_user)
    end if

    if (present(interp_type2)) then
      select case (trim(interp_type2))
        case ('line', 'linear')
          call cvmix_put_kpp('interp_type2', CVMIX_MATH_INTERP_LINEAR,        &
                             CVmix_kpp_params_user)
        case ('quad', 'quadratic')
          call cvmix_put_kpp('interp_type2', CVMIX_MATH_INTERP_QUAD,          &
                             CVmix_kpp_params_user)
        case ('cube', 'cubic', 'cubic_spline', 'cubic spline')
          call cvmix_put_kpp('interp_type2', CVMIX_MATH_INTERP_CUBE_SPLINE,   &
                             CVmix_kpp_params_user)
        case ('POP','LMD94')
          call cvmix_put_kpp('interp_type2', CVMIX_KPP_INTERP_LMD94,          &
                             CVmix_kpp_params_user)
        case DEFAULT
          print*, "ERROR: ", trim(interp_type2), " is not a valid type of ",  &
                  "interpolation!"
          stop 1
      end select
    else
      call cvmix_put_kpp('interp_type2', CVMIX_KPP_INTERP_LMD94,              &
                         CVmix_kpp_params_user)
    end if

    if (present(MatchTechnique)) then
      select case (trim(MatchTechnique))
        case ('MatchBoth')
          call cvmix_put_kpp('MatchTechnique', CVMIX_KPP_MATCH_BOTH,          &
                             CVmix_kpp_params_user)
        case ('MatchGradient')
          call cvmix_put_kpp('MatchTechnique', CVMIX_KPP_MATCH_GRADIENT,      &
                             CVmix_kpp_params_user)
        case ('SimpleShapes')
          call cvmix_put_kpp('MatchTechnique', CVMIX_KPP_SIMPLE_SHAPES,       &
                             CVmix_kpp_params_user)
        case ('ParabolicNonLocal')
          call cvmix_put_kpp('MatchTechnique', CVMIX_KPP_PARABOLIC_NONLOCAL,  &
                             CVmix_kpp_params_user)
        case DEFAULT
          print*, "ERROR: ", trim(MatchTechnique), " is not a valid choice ", &
                  "for MatchTechnique!"
          stop 1
        end select
    else
      call cvmix_put_kpp('MatchTechnique', CVMIX_KPP_SIMPLE_SHAPES,           &
                         CVmix_kpp_params_user)
    end if

    if (present(old_vals)) then
      select case (trim(old_vals))
        case ("overwrite")
          call cvmix_put_kpp('handle_old_vals', CVMIX_OVERWRITE_OLD_VAL,      &
                               cvmix_kpp_params_user)
        case ("sum")
          call cvmix_put_kpp('handle_old_vals', CVMIX_SUM_OLD_AND_NEW_VALS,   &
                               cvmix_kpp_params_user)
        case ("max")
          call cvmix_put_kpp('handle_old_vals', CVMIX_MAX_OLD_AND_NEW_VALS,   &
                               cvmix_kpp_params_user)
        case DEFAULT
          print*, "ERROR: ", trim(old_vals), " is not a valid option for ",   &
                  "handling old values of diff and visc."
          stop 1
      end select
    else
      call cvmix_put_kpp('handle_old_vals', CVMIX_OVERWRITE_OLD_VAL,          &
                               cvmix_kpp_params_user)
    end if

    if (present(lStokesMOST)) then
       call cvmix_put_kpp('lStokesMOST', lStokesMOST, CVmix_kpp_params_user)
    else
      call cvmix_put_kpp('lStokesMOST', .false., CVmix_kpp_params_user)
    end if
    call cvmix_put_kpp('CVt2', 2.6_cvmix_r8, CVmix_kpp_params_user)

    if (present(lEkman)) then
      call cvmix_put_kpp('lEkman', lEkman, CVmix_kpp_params_user)
    else
      call cvmix_put_kpp('lEkman', .false., CVmix_kpp_params_user)
    end if

    if (present(lMonOb)) then
      call cvmix_put_kpp('lMonOb', lMonOb, CVmix_kpp_params_user)
    else
      call cvmix_put_kpp('lMonOb', .false., CVmix_kpp_params_user)
    end if

    if (present(lnoDGat1)) then
      call cvmix_put_kpp('lnoDGat1', lnoDGat1, CVmix_kpp_params_user)
    else
      call cvmix_put_kpp('lnoDGat1', .true., CVmix_kpp_params_user)
    end if

    if (present(lenhanced_diff)) then
      call cvmix_put_kpp('lenhanced_diff', lenhanced_diff,                    &
                         CVmix_kpp_params_user)
    else
      call cvmix_put_kpp('lenhanced_diff', .true., CVmix_kpp_params_user)
    end if

    if (present(Langmuir_mixing_str)) then
       select case (trim(Langmuir_mixing_str))
       case ("LWF16")
          call cvmix_put_kpp('Langmuir_Mixing_Opt', LANGMUIR_MIXING_LWF16 ,   &
               CVmix_kpp_params_user)
       case ("RWHGK16")
          call cvmix_put_kpp('Langmuir_Mixing_Opt',                           &
               LANGMUIR_MIXING_RWHGK16, CVmix_kpp_params_user)
       case ("NONE")
          call cvmix_put_kpp('Langmuir_Mixing_Opt',                           &
               NO_LANGMUIR_MIXING, CVmix_kpp_params_user)
       case DEFAULT
          print*, "ERROR: ", trim(Langmuir_mixing_str), " is not a valid ",   &
                  "option for Langmuir_mixing_str!"
          stop 1
       end select
    else
       call cvmix_put_kpp('Langmuir_Mixing_Opt',                              &
            NO_LANGMUIR_MIXING, CVmix_kpp_params_user)
    end if

    if (present(Langmuir_entrainment_str)) then
       select case (trim(Langmuir_entrainment_str))
       case ("LWF16")
          call cvmix_put_kpp('Langmuir_Entrainment_Opt',                      &
               LANGMUIR_ENTRAINMENT_LWF16, CVmix_kpp_params_user)
       case ("LF17")
          call cvmix_put_kpp('Langmuir_Entrainment_Opt',                      &
               LANGMUIR_ENTRAINMENT_LF17, CVmix_kpp_params_user)
       case ("RWHGK16")
          call cvmix_put_kpp('Langmuir_Entrainment_Opt',                      &
               LANGMUIR_ENTRAINMENT_RWHGK16, CVmix_kpp_params_user)
       case ("NONE")
          call cvmix_put_kpp('Langmuir_Entrainment_Opt',                      &
               NO_LANGMUIR_ENTRAINMENT, CVmix_kpp_params_user)
       case DEFAULT
          print*, "ERROR: ", trim(Langmuir_entrainment_str), " is not a ",    &
                  "valid option for Langmuir_entrainment_str!"
          stop 1
       end select
    else
       call cvmix_put_kpp('Langmuir_Entrainment_Opt',                         &
            NO_LANGMUIR_ENTRAINMENT, CVmix_kpp_params_user)
    end if

    ! By default, assume that G(0) = 0 for nonlocal term
    nonlocal_coeff = (Cstar_loc*vonkar_loc*                                   &
                      (vonkar_loc*surf_layer_ext_loc*c_s)**                   &
                      (cvmix_one/real(3,cvmix_r8)))
    if (present(lnonzero_surf_nonlocal)) then
      if (lnonzero_surf_nonlocal) then
        nonlocal_coeff = real(1,cvmix_r8)
      end if
    end if
    call cvmix_put_kpp('nonlocal_coeff',nonlocal_coeff,CVmix_kpp_params_user)

    ! By default, use sigma construction from Danabasoglu et al. when computing
    ! turbulent scales. Set l_LMD_ws = .true. to use Large et al. construction.
    if (present(l_LMD_ws)) then
      call cvmix_put_kpp('l_LMD_ws', l_LMD_ws,                    &
                         CVmix_kpp_params_user)
    else
      call cvmix_put_kpp('l_LMD_ws', .false., CVmix_kpp_params_user)
    end if

    ! Initialize parameters for enhanced entrainment
    call cvmix_put_kpp('c_ST', 0.17_cvmix_r8, CVmix_kpp_params_user)
    call cvmix_put_kpp('c_CT', 0.15_cvmix_r8, CVmix_kpp_params_user)
    call cvmix_put_kpp('c_LT', 0.083_cvmix_r8, CVmix_kpp_params_user)
    call cvmix_put_kpp('p_LT', 2.0_cvmix_r8, CVmix_kpp_params_user)
    call cvmix_put_kpp('RWHGK_ENTR_COEF', 2.3_cvmix_r8, CVmix_kpp_params_user)
    call cvmix_put_kpp('RWHGK_ENTR_EXP', -0.5_cvmix_r8, CVmix_kpp_params_user)
!EOC

  end subroutine cvmix_init_kpp

!BOP

! !IROUTINE: cvmix_coeffs_kpp_wrap
! !INTERFACE:

  subroutine cvmix_coeffs_kpp_wrap(CVmix_vars, CVmix_kpp_params_user)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for the KPP boundary layer mixing
!  parameterization.
!\\
!\\

! !INPUT PARAMETERS:
    type(cvmix_kpp_params_type), intent(in), optional, target ::              &
                                           CVmix_kpp_params_user

! !INPUT/OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars

!EOP
!BOC

    real(cvmix_r8), dimension(CVmix_vars%max_nlev+1) :: new_Mdiff, new_Tdiff, &
                                                        new_Sdiff
    integer :: nlev, max_nlev
    type(cvmix_kpp_params_type),  pointer :: CVmix_kpp_params_in

    CVmix_kpp_params_in => CVmix_kpp_params_saved
    if (present(CVmix_kpp_params_user)) then
      CVmix_kpp_params_in => CVmix_kpp_params_user
    end if

    nlev = CVmix_vars%nlev
    max_nlev = CVmix_vars%max_nlev

    if (.not.associated(CVmix_vars%Mdiff_iface)) &
      call cvmix_put(CVmix_vars, "Mdiff", cvmix_zero, max_nlev)
    if (.not.associated(CVmix_vars%Tdiff_iface)) &
      call cvmix_put(CVmix_vars, "Tdiff", cvmix_zero, max_nlev)
    if (.not.associated(CVmix_vars%Sdiff_iface)) &
      call cvmix_put(CVmix_vars, "Sdiff", cvmix_zero, max_nlev)

    call cvmix_put(CVmix_vars, 'kpp_transport', cvmix_zero, max_nlev)

    call cvmix_coeffs_kpp(new_Mdiff, new_Tdiff, new_Sdiff,                    &
                          CVmix_vars%zw_iface, CVmix_vars%zt_cntr,            &
                          CVmix_vars%Mdiff_iface, CVmix_vars%Tdiff_iface,     &
                          CVMix_vars%Sdiff_iface,                             &
                          CVmix_vars%BoundaryLayerDepth,                      &
                          CVmix_vars%kOBL_depth,                              &
                          CVmix_vars%kpp_Tnonlocal_iface,                     &
                          CVmix_vars%kpp_Snonlocal_iface,                     &
                          CVmix_vars%SurfaceFriction,                         &
                          CVmix_vars%SurfaceBuoyancyForcing,                  &
                          nlev, max_nlev,                                     &
                          CVmix_vars%LangmuirEnhancementFactor,               &
                          CVmix_vars%StokesMostXi,                            &
                          CVmix_kpp_params_user)

    call cvmix_update_wrap(CVmix_kpp_params_in%handle_old_vals, max_nlev,     &
                           Mdiff_out = CVmix_vars%Mdiff_iface,                &
                           new_Mdiff = new_Mdiff,                             &
                           Tdiff_out = CVmix_vars%Tdiff_iface,                &
                           new_Tdiff = new_Tdiff,                             &
                           Sdiff_out = CVmix_vars%Sdiff_iface,                &
                           new_Sdiff = new_Sdiff)

!EOC

  end subroutine cvmix_coeffs_kpp_wrap

!BOP

! !IROUTINE: cvmix_coeffs_kpp_low
! !INTERFACE:

  subroutine cvmix_coeffs_kpp_low(Mdiff_out, Tdiff_out, Sdiff_out, zw, zt,    &
                                  old_Mdiff, old_Tdiff, old_Sdiff, OBL_depth, &
                                  kOBL_depth, Tnonlocal, Snonlocal, surf_fric,&
                                  surf_buoy, nlev, max_nlev, Langmuir_EFactor,&
                                  StokesXI,CVmix_kpp_params_user)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for the KPP boundary layer mixing
!  parameterization.
!\\
!\\

! !INPUT PARAMETERS:
    type(cvmix_kpp_params_type),  intent(in), optional, target ::             &
                                            CVmix_kpp_params_user

    integer,                               intent(in) :: nlev, max_nlev
    real(cvmix_r8), dimension(max_nlev+1), intent(in) :: old_Mdiff,           &
                                                         old_Tdiff,           &
                                                         old_Sdiff,           &
                                                         zw
    real(cvmix_r8), dimension(max_nlev),   intent(in) :: zt
    real(cvmix_r8),                        intent(in) :: OBL_depth,           &
                                                         surf_fric,           &
                                                         surf_buoy,           &
                                                         kOBL_depth
    ! Langmuir enhancement factor
    real(cvmix_r8), intent(in), optional :: Langmuir_EFactor
    real(cvmix_r8), intent(in), optional :: StokesXI
! !INPUT/OUTPUT PARAMETERS:
    real(cvmix_r8), dimension(max_nlev+1), intent(inout) :: Mdiff_out,        &
                                                            Tdiff_out,        &
                                                            Sdiff_out,        &
                                                            Tnonlocal,        &
                                                            Snonlocal

!EOP
!BOC

    ! Local variables
    type(cvmix_kpp_params_type), pointer :: CVmix_kpp_params_in

    ! OBL_[MTS]diff are the diffusivities in the whole OBL
    real(cvmix_r8), dimension(nint(kOBL_depth)) :: OBL_Mdiff, OBL_Tdiff,      &
                                                   OBL_Sdiff

    ! [MTS]diff_ktup are the enhanced diffusivity and viscosity values at the
    ! deepest cell center above OBL_depth. Other _ktup vars are intermediary
    ! variables needed to compute [MTS]diff_ktup
    real(cvmix_r8) :: Mdiff_ktup, Tdiff_ktup, Sdiff_ktup
    real(cvmix_r8) :: sigma_ktup, wm_ktup, ws_ktup

    real(cvmix_r8) :: delta

    real(cvmix_r8), dimension(nlev+1) :: sigma, w_m, w_s

    ! [MTS]shape are the coefficients of the shape function in the gradient
    ! term; [TS]shape2 are the coefficients for the nonlocal term;
    ! NMshape is the coefficient for the no-matching case, an option to shape
    ! a Langmuir enhancement
    real(cvmix_r8), dimension(4) :: Mshape, Tshape, Sshape, Tshape2, Sshape2,&
                                    NMshape

    ! [MTS]shapeAt1 is value of shape function at sigma = 1
    ! d[MTS]shapeAt1 is value of derivative of shape function at sigma = 1
    ! (Used for matching the shape function at OBL depth)
    real(cvmix_r8) :: MshapeAt1, TshapeAt1, SshapeAt1
    real(cvmix_r8) :: dMshapeAt1, dTshapeAt1, dSshapeAt1

    ! [MTS]shapeAtS is value of shape function at sigma = S
    real(cvmix_r8) :: MshapeAtS, TshapeAtS, SshapeAtS, GAtS
    ! Storing the maximum value of shape function for no-matching case
    !  that is used as an option for Langmuir mixing
    real(cvmix_r8), parameter :: NMshapeMax = 4./27.

    ! [MTS]diff_OBL is value of diffusivity at OBL depth
    ! d[MTS]diff_OBL is value of derivative of diffusivity at OBL depth
    ! w[ms]_OBL is value of wm or ws at OBL depth
    real(cvmix_r8) :: Mdiff_OBL, Tdiff_OBL, Sdiff_OBL
    real(cvmix_r8) :: dMdiff_OBL, dTdiff_OBL, dSdiff_OBL
    real(cvmix_r8) :: wm_OBL, ws_OBL, second_term

    ! coefficients used for interpolation if interp_type2 is not 'LMD94'
    real(kind=cvmix_r8), dimension(4) :: coeffs

    ! Width of column kw_up and kw_up+1
    real(cvmix_r8), dimension(2) :: col_widths, col_centers
    real(cvmix_r8), dimension(2) :: Mdiff_vals, Tdiff_vals, Sdiff_vals

    ! Parameters for RWHGK16 Langmuir parameterization
    real(cvmix_r8) :: MixingCoefEnhancement
    real(cvmix_r8) :: ShapeNoMatchAtS

    ! Parameters for Stokes_MOST
    real(cvmix_r8) :: Gcomposite, Hsigma, sigh, T_NLenhance , S_NLenhance , XIone

    ! Constant from params
    integer :: interp_type2, MatchTechnique

    integer :: kw
    logical :: lstable
    integer :: ktup, & ! kt index of cell center above OBL_depth
               kwup    ! kw index of iface above OBL_depth (= kt index of
                       ! cell containing OBL_depth)

    CVmix_kpp_params_in => CVmix_kpp_params_saved
    if (present(CVmix_kpp_params_user)) then
      CVmix_kpp_params_in => CVmix_kpp_params_user
    end if
    interp_type2   = CVmix_kpp_params_in%interp_type2
    MatchTechnique = CVmix_kpp_params_in%MatchTechnique

    ! Output values should be set to input values
    Mdiff_out = old_Mdiff
    Tdiff_out = old_Tdiff
    Sdiff_out = old_Sdiff

    ! (1) Column-specific parameters
    !
    ! Stability => positive surface buoyancy flux
    lstable = (surf_buoy.gt.cvmix_zero)

    kwup = floor(kOBL_depth)
    ktup = nint(kOBL_depth)-1

    if (ktup.eq.nlev) then
      ! OBL_depth between bottom cell center and ocean bottom, assume
      ! zt(ktup+1) = ocn_bottom (which is zw(nlev+1)
      delta = (OBL_depth+zt(ktup))/(zt(ktup)-zw(ktup+1))
    else
      delta = (OBL_depth+zt(ktup))/(zt(ktup)-zt(ktup+1))
    end if

  if ( CVmix_kpp_params_in%lStokesMOST ) then             ! Stokes_MOST

        ! (2a) Compute turbulent scales at OBL depth
        call cvmix_kpp_compute_turbulent_scales(cvmix_one, OBL_depth,         &
                                                surf_buoy, surf_fric,         &
                               StokesXI,        wm_OBL, ws_OBL,               &
                                                CVmix_kpp_params_user)

        ! (2b) Compute diffusivities at OBL depth
          col_centers(1) = zt(kwup)
          col_widths(1) = zw(kwup) - zw(kwup+1)
          Mdiff_vals(1) = old_Mdiff(kwup+1)
          Tdiff_vals(1) = old_Tdiff(kwup+1)
          Sdiff_vals(1) = old_Sdiff(kwup+1)
          if (kwup.eq.nlev) then
            col_centers(2) = zw(kwup+1)
            col_widths(2)  = 1.0_cvmix_r8 !Value doesn't matter, will divide into zero
            Mdiff_vals(2)  = old_Mdiff(kwup+1)   !Mdiff_out(kwup+1)
            Tdiff_vals(2)  = old_Tdiff(kwup+1)
            Sdiff_vals(2)  = old_Sdiff(kwup+1)
          else
            col_centers(2) = zt(kwup+1)
            col_widths(2)  = zw(kwup+1) - zw(kwup+2)
            Mdiff_vals(2)  = old_Mdiff(kwup+2)
            Tdiff_vals(2)  = old_Tdiff(kwup+2)
            Sdiff_vals(2)  = old_Sdiff(kwup+2)
          end if
          if (kwup.eq.1) then
            Mdiff_OBL = cvmix_kpp_compute_nu_at_OBL_depth_LMD94(col_centers,  &
                                                      col_widths,             &
                                                      Mdiff_vals, OBL_depth,  &
                                                      dnu_dz=dMdiff_OBL)
            Tdiff_OBL = cvmix_kpp_compute_nu_at_OBL_depth_LMD94(col_centers,  &
                                                      col_widths,             &
                                                      Tdiff_vals, OBL_depth,  &
                                                      dnu_dz=dTdiff_OBL)
            Sdiff_OBL = cvmix_kpp_compute_nu_at_OBL_depth_LMD94(col_centers,  &
                                                      col_widths,             &
                                                      Sdiff_vals, OBL_depth,  &
                                                      dnu_dz=dSdiff_OBL)
          else ! interp_type == 'LMD94' and kwup > 1
            Mdiff_OBL = cvmix_kpp_compute_nu_at_OBL_depth_LMD94(col_centers,  &
                                                      col_widths,             &
                                                      Mdiff_vals, OBL_depth,  &
                                                      old_Mdiff(kwup),        &
                                                      dnu_dz=dMdiff_OBL)
            Tdiff_OBL = cvmix_kpp_compute_nu_at_OBL_depth_LMD94(col_centers,  &
                                                      col_widths,             &
                                                      Tdiff_vals, OBL_depth,  &
                                                      old_Tdiff(kwup),        &
                                                      dnu_dz=dTdiff_OBL)
            Sdiff_OBL = cvmix_kpp_compute_nu_at_OBL_depth_LMD94(col_centers,  &
                                                      col_widths,             &
                                                      Sdiff_vals, OBL_depth,  &
                                                      old_Sdiff(kwup),        &
                                                      dnu_dz=dSdiff_OBL)
          end if

    !     (3a) Compute turbulent scales at interfaces throughout column
    sigma = -zw(1:nlev+1)/OBL_depth
    call cvmix_kpp_compute_turbulent_scales(sigma, OBL_depth, surf_buoy, & !_1d
                surf_fric,   xi=StokesXI,    w_m=w_m, w_s=w_s,           &
                CVmix_kpp_params_user = CVmix_kpp_params_user)

    do kw=2,kwup                                   ! OBL overwrite loop to kwup
    !     (3b) Evaluate G(sigma) >= 0  at each cell interface
      Gcomposite = cvmix_kpp_composite_shape(sigma(kw))
      sigh   = MAX(CVmix_kpp_params_in%surf_layer_ext, MIN(sigma(kw) ,cvmix_one))
      Hsigma = ((sigh      - CVmix_kpp_params_in%surf_layer_ext) / &
                (cvmix_one - CVmix_kpp_params_in%surf_layer_ext) )**2
!     Hsigma = MAX( cvmix_zero , MIN( cvmix_one , Hsigma ) )

    !   (3c) Compute nonlocal term at each cell interface
      if (.not.lstable) then
        Tnonlocal(kw) = 4.7 * Gcomposite ! LMD 6.26 Gcubic
        Snonlocal(kw) = 4.7 * Gcomposite
      else
        Tnonlocal(kw) = cvmix_zero
        Snonlocal(kw) = cvmix_zero
      end if

    !   (3d) Diffusivity = (OBL_depth * turbulent scale * G(sigma) + Xdiff_OBL * Hsigma)
      OBL_Mdiff(kw) = OBL_depth * w_m(kw) * Gcomposite  + Mdiff_OBL * Hsigma
      OBL_Tdiff(kw) = OBL_depth * w_s(kw) * Gcomposite  + Tdiff_OBL * Hsigma
      OBL_Sdiff(kw) = OBL_depth * w_s(kw) * Gcomposite  + Sdiff_OBL * Hsigma

    end do

    ! (4) Compute the enhanced diffusivity
    !     (4a) Compute shape function at last cell center in OBL
    sigma_ktup = -zt(ktup)/OBL_depth
    Gcomposite =  cvmix_kpp_composite_shape(sigma_ktup)
    sigh   = MAX(CVmix_kpp_params_in%surf_layer_ext, MIN(sigma_ktup ,cvmix_one))
    Hsigma = ( (sigh - CVmix_kpp_params_in%surf_layer_ext) / &
                 (cvmix_one - CVmix_kpp_params_in%surf_layer_ext) )**2

    !     (4b) Compute turbulent scales at last cell center in OBL
    call cvmix_kpp_compute_turbulent_scales(sigma_ktup, OBL_depth, surf_buoy, & !0d
                      surf_fric, StokesXI,  wm_ktup, ws_ktup,      &
                     CVmix_kpp_params_user)

    !     (4c) Diffusivity at last cell center in OBL
    Mdiff_ktup = OBL_depth * wm_ktup * Gcomposite  + Mdiff_OBL * Hsigma
    Tdiff_ktup = OBL_depth * ws_ktup * Gcomposite  + Tdiff_OBL * Hsigma
    Sdiff_ktup = OBL_depth * ws_ktup * Gcomposite  + Sdiff_OBL * Hsigma

    if (CVmix_kpp_params_in%lenhanced_diff) then
      if ((ktup.eq.kwup).or.(ktup.eq.kwup-1)) then
        T_NLenhance =  Tnonlocal(ktup+1)
        S_NLenhance =  Snonlocal(ktup+1)
        call cvmix_kpp_compute_enhanced_diff(Mdiff_ktup,                      &
                                             Tdiff_ktup,                      &
                                             Sdiff_ktup,                      &
                                             Mdiff_out(ktup+1),               &
                                             Tdiff_out(ktup+1),               &
                                             Sdiff_out(ktup+1),               &
                                             OBL_Mdiff(ktup+1),               &
                                             OBL_Tdiff(ktup+1),               &
                                             OBL_Sdiff(ktup+1),               &
                                             T_NLenhance      ,               &
                                             S_NLenhance      ,               &
                                             delta, lkteqkw=(ktup.eq.kwup))
      else
        print*, "ERROR: ktup should be either kwup or kwup-1!"
        print*, "ktup = ", ktup, " and kwup = ", kwup
        stop 1
      end if
    else
      if ( kwup .eq. ktup )  then
        OBL_Mdiff(ktup+1) = old_Mdiff(ktup+1)
        OBL_Tdiff(ktup+1) = old_Tdiff(ktup+1)
        OBL_Sdiff(ktup+1) = old_Sdiff(ktup+1)
      end if
    end if

    ! (5) Combine interior and boundary coefficients

    Mdiff_out(2:ktup+1) = OBL_Mdiff(2:ktup+1)
    Tdiff_out(2:ktup+1) = OBL_Tdiff(2:ktup+1)
    Sdiff_out(2:ktup+1) = OBL_Sdiff(2:ktup+1)
  else ! not  Stokes_MOST

    XIone = cvmix_one
    ! (2) Compute coefficients of shape function
    !     A no-match case is stored for use in Langmuir scheme
    NMshape(1) = cvmix_zero
    NMshape(2) = cvmix_one
    NMshape(3) = -real(2,cvmix_r8)
    NMshape(4) = cvmix_one
    select case (MatchTechnique)
      case (CVMIX_KPP_SIMPLE_SHAPES)
        ! Simple shape function is sigma*(1-sigma)^2
        Mshape(1) =  cvmix_zero
        Mshape(2) =  cvmix_one
        Mshape(3) = -real(2,cvmix_r8)
        Mshape(4) =  cvmix_one
        Tshape    = Mshape
        Sshape    = Mshape
        Tshape2   = Tshape
        Sshape2   = Sshape
      case (CVMIX_KPP_PARABOLIC_NONLOCAL)
        ! Shape function is sigma*(1-sigma)^2 for gradient term
        ! and (1-sigma)^2 for non-local term
        Mshape(1) =  cvmix_zero
        Mshape(2) =  cvmix_one
        Mshape(3) = -real(2,cvmix_r8)
        Mshape(4) =  cvmix_one
        Tshape    = Mshape
        Sshape    = Mshape
        Tshape2(1) =  cvmix_one
        Tshape2(2) = -real(2,cvmix_r8)
        Tshape2(3) =  cvmix_one
        Tshape2(4) =  cvmix_zero
        Sshape2    = Tshape2
      case DEFAULT
        ! (2a) Compute turbulent scales at OBL depth
        call cvmix_kpp_compute_turbulent_scales(cvmix_one, OBL_depth,         &
                                                surf_buoy, surf_fric,         &
                                                XIone, wm_OBL, ws_OBL,        &
                                                CVmix_kpp_params_user)
        if (CVMix_KPP_Params_in%Langmuir_Mixing_Opt &
           .eq. LANGMUIR_MIXING_LWF16) then
          ! enhance the turbulent velocity scale
          wm_OBL = wm_OBL * Langmuir_EFactor
          ws_OBL = ws_OBL * Langmuir_EFactor
        end if
        ! (2b) Compute diffusivities at OBL depth
        if (interp_type2.ne.CVMIX_KPP_INTERP_LMD94) then
          if (kwup.eq.1) then
            call cvmix_math_poly_interp(coeffs, interp_type2, zw(kwup:kwup+1),&
                                        old_Mdiff(kwup:kwup+1))
            Mdiff_OBL = cvmix_math_evaluate_cubic(coeffs, -OBL_depth,         &
                                                  dMdiff_OBL)

            call cvmix_math_poly_interp(coeffs, interp_type2, zw(kwup:kwup+1),&
                                        old_Tdiff(kwup:kwup+1))
            Tdiff_OBL = cvmix_math_evaluate_cubic(coeffs, -OBL_depth,         &
                                                  dTdiff_OBL)

            call cvmix_math_poly_interp(coeffs, interp_type2, zw(kwup:kwup+1),&
                                        old_Sdiff(kwup:kwup+1))
            Sdiff_OBL = cvmix_math_evaluate_cubic(coeffs, -OBL_depth,         &
                                                  dSdiff_OBL)
          else ! interp_type2 != 'LMD94' and kwup > 1
            call cvmix_math_poly_interp(coeffs, interp_type2, zw(kwup:kwup+1),&
                                        old_Mdiff(kwup:kwup+1), zw(kwup-1),   &
                                        old_Mdiff(kwup-1))
            Mdiff_OBL = cvmix_math_evaluate_cubic(coeffs, -OBL_depth,         &
                                                  dMdiff_OBL)

            call cvmix_math_poly_interp(coeffs, interp_type2, zw(kwup:kwup+1),&
                                        old_Tdiff(kwup:kwup+1), zw(kwup-1),   &
                                        old_Tdiff(kwup-1))
            Tdiff_OBL = cvmix_math_evaluate_cubic(coeffs, -OBL_depth,         &
                                                  dTdiff_OBL)

            call cvmix_math_poly_interp(coeffs, interp_type2, zw(kwup:kwup+1),&
                                        old_Sdiff(kwup:kwup+1), zw(kwup-1),   &
                                        old_Sdiff(kwup-1))
            Sdiff_OBL = cvmix_math_evaluate_cubic(coeffs, -OBL_depth,         &
                                                  dSdiff_OBL)
          end if
        else ! interp_type2 == 'LMD94'
          col_centers(1) = zt(kwup)
          col_widths(1) = zw(kwup) - zw(kwup+1)
          Mdiff_vals(1) = old_Mdiff(kwup+1)
          Tdiff_vals(1) = old_Tdiff(kwup+1)
          Sdiff_vals(1) = old_Sdiff(kwup+1)
          if (kwup.eq.nlev) then
            col_centers(2) = zw(kwup+1)
            col_widths(2)  = 1.0_cvmix_r8 ! Value doesn't matter, will divide
                                          ! into zero
            Mdiff_vals(2)  = old_Mdiff(kwup+1)
            Tdiff_vals(2)  = old_Tdiff(kwup+1)
            Sdiff_vals(2)  = old_Sdiff(kwup+1)
          else
            col_centers(2) = zt(kwup+1)
            col_widths(2)  = zw(kwup+1) - zw(kwup+2)
            Mdiff_vals(2)  = old_Mdiff(kwup+2)
            Tdiff_vals(2)  = old_Tdiff(kwup+2)
            Sdiff_vals(2)  = old_Sdiff(kwup+2)
          end if

          if (kwup.eq.1) then
            Mdiff_OBL = cvmix_kpp_compute_nu_at_OBL_depth_LMD94(col_centers,  &
                                                      col_widths,             &
                                                      Mdiff_vals, OBL_depth,  &
                                                      dnu_dz=dMdiff_OBL)
            Tdiff_OBL = cvmix_kpp_compute_nu_at_OBL_depth_LMD94(col_centers,  &
                                                      col_widths,             &
                                                      Tdiff_vals, OBL_depth,  &
                                                      dnu_dz=dTdiff_OBL)
            Sdiff_OBL = cvmix_kpp_compute_nu_at_OBL_depth_LMD94(col_centers,  &
                                                      col_widths,             &
                                                      Sdiff_vals, OBL_depth,  &
                                                      dnu_dz=dSdiff_OBL)
          else ! interp_type == 'LMD94' and kwup > 1
            Mdiff_OBL = cvmix_kpp_compute_nu_at_OBL_depth_LMD94(col_centers,  &
                                                      col_widths,             &
                                                      Mdiff_vals, OBL_depth,  &
                                                      old_Mdiff(kwup),        &
                                                      dnu_dz=dMdiff_OBL)
            Tdiff_OBL = cvmix_kpp_compute_nu_at_OBL_depth_LMD94(col_centers,  &
                                                      col_widths,             &
                                                      Tdiff_vals, OBL_depth,  &
                                                      old_Tdiff(kwup),        &
                                                      dnu_dz=dTdiff_OBL)
            Sdiff_OBL = cvmix_kpp_compute_nu_at_OBL_depth_LMD94(col_centers,  &
                                                      col_widths,             &
                                                      Sdiff_vals, OBL_depth,  &
                                                      old_Sdiff(kwup),        &
                                                      dnu_dz=dSdiff_OBL)
          end if
        end if ! interp_type != "LMD94"

        ! (2c) Compute G(1) [shape function when sigma = 1] and G'(1) for three
        !      cases:
        !      i) momentum diffusivity (viscosity)
        !      ii) temperature diffusivity
        !      iii) other tracers diffusivity
        ! Notes:
        !   * We are computing G(1) and G'(1) so we can represent G(sigma) as a
        !     cubic polynomial and then compute Kx = OBL_depth*wx*G. If either
        !     OBL_depth or wx are 0, it doesn't matter what G is because Kx
        !     will be zero everywhere... in these cases, we set G(1)=G'(1)=0.
        !   * If OBL_depth = 0, the above note applies to all three situations
        !     listed as (i), (ii), and (iii). If ws = 0, it applies only to (i)
        !     and (ii). If wm = 0, it applies only to (iii).
        if (OBL_depth.eq.cvmix_zero) then
          ! Values don't matter, K = 0
          MshapeAt1 = cvmix_zero
          TshapeAt1 = cvmix_zero
          SshapeAt1 = cvmix_zero
          dMshapeAt1 = cvmix_zero
          dTshapeAt1 = cvmix_zero
          dSshapeAt1 = cvmix_zero
        else ! OBL_depth != 0
          if (wm_OBL.ne.cvmix_zero) then
            MshapeAt1 = Mdiff_OBL/(wm_OBL*OBL_depth)
          else
            MshapeAt1 = cvmix_zero ! value doesn't really matter, Km = 0
          end if
          if (ws_OBL.ne.cvmix_zero) then
            TshapeAt1 = Tdiff_OBL/(ws_OBL*OBL_depth)
            SshapeAt1 = Sdiff_OBL/(ws_OBL*OBL_depth)
          else
            TshapeAt1 = cvmix_zero ! value doesn't really matter, Ks = 0
            SshapeAt1 = cvmix_zero ! value doesn't really matter, Ks = 0
          end if
          if (CVmix_kpp_params_in%lnoDGat1) then
            ! Force G'(1) = 0
            dMshapeAt1 = cvmix_zero
            dTshapeAt1 = cvmix_zero
            dSshapeAt1 = cvmix_zero
          else
            second_term = real(5,cvmix_r8)*surf_buoy/(surf_fric**4)
            if (wm_OBL.ne.cvmix_zero) then
              dMshapeAt1 = -dMdiff_OBL/wm_OBL
              if (lstable) &
                dMshapeAt1 = dMshapeAt1 + second_term*Mdiff_OBL
            else
              dMshapeAt1 = cvmix_zero ! value doesn't really matter, Km = 0
            end if
            if (ws_OBL.ne.cvmix_zero) then
              dTshapeAt1 = -dTdiff_OBL/ws_OBL
              dSshapeAt1 = -dSdiff_OBL/ws_OBL
              if (lstable) then
                dTshapeAt1 = dTshapeAt1 + second_term*Tdiff_OBL
                dSshapeAt1 = dSshapeAt1 + second_term*Sdiff_OBL
              end if
            else
              dTshapeAt1 = cvmix_zero ! value doesn't really matter, Ks = 0
              dSshapeAt1 = cvmix_zero ! value doesn't really matter, Ks = 0
            end if
            dMshapeAt1 = min(dMshapeAt1, cvmix_zero) ! non-positive value!
            dTshapeAt1 = min(dTshapeAt1, cvmix_zero) ! non-positive value!
            dSshapeAt1 = min(dSshapeAt1, cvmix_zero) ! non-positive value!
          end if ! lnoDGat1
        end if ! OBL_depth == 0

        !   (2d) Compute coefficients of shape function
        call cvmix_kpp_compute_shape_function_coeffs(MshapeAt1, dMshapeAt1,   &
                                                     Mshape)
        call cvmix_kpp_compute_shape_function_coeffs(TshapeAt1, dTshapeAt1,   &
                                                     Tshape)
        call cvmix_kpp_compute_shape_function_coeffs(SshapeAt1, dSshapeAt1,   &
                                                     Sshape)
        if (MatchTechnique.eq.CVMIX_KPP_MATCH_GRADIENT) then
          ! Only match for gradient term, use simple shape for nonlocal
          Tshape2(1) =  cvmix_zero
          Tshape2(2) =  cvmix_one
          Tshape2(3) = -real(2,cvmix_r8)
          Tshape2(4) =  cvmix_one
          Sshape2 = Tshape2
        else
          ! Shape function is the same for gradient and nonlocal
          Tshape2 = Tshape
          Sshape2 = Sshape
        end if
    end select

    ! (3) Use shape function to compute diffusivities throughout OBL
    Tnonlocal = cvmix_zero
    Snonlocal = cvmix_zero
    OBL_Mdiff = cvmix_zero
    OBL_Tdiff = cvmix_zero
    OBL_Sdiff = cvmix_zero
    sigma = -zw(1:nlev+1)/OBL_depth
    !     (3a) Compute turbulent scales throghout column
    call cvmix_kpp_compute_turbulent_scales(sigma, OBL_depth, surf_buoy,   &
                                            surf_fric, XIone, w_m, w_s,    &
                                            CVmix_kpp_params_user)
    do kw=2,kwup
      !   (3b) Evaluate G(sigma) at each cell interface
      MshapeAtS = cvmix_math_evaluate_cubic(Mshape, sigma(kw))
      TshapeAtS = cvmix_math_evaluate_cubic(Tshape, sigma(kw))
      SshapeAtS = cvmix_math_evaluate_cubic(Sshape, sigma(kw))
      ! The RWHGK16 Langmuir uses the shape function to shape the
      !  enhancement to the mixing coefficient.
      ShapeNoMatchAtS = cvmix_math_evaluate_cubic(NMshape, sigma(kw))
      !   (3c) Compute nonlocal term at each cell interface
      if (.not.lstable) then
        GAtS = cvmix_math_evaluate_cubic(Tshape2, sigma(kw))
        Tnonlocal(kw) = CVmix_kpp_params_in%nonlocal_coeff*GAtS
        GAtS = cvmix_math_evaluate_cubic(Sshape2, sigma(kw))
        Snonlocal(kw) = CVmix_kpp_params_in%nonlocal_coeff*GAtS
      end if

      select case (CVMix_KPP_Params_in%Langmuir_Mixing_Opt)
      case (LANGMUIR_MIXING_LWF16)
         MixingCoefEnhancement = Langmuir_EFactor
      case (LANGMUIR_MIXING_RWHGK16)
         MixingCoefEnhancement = cvmix_one + ShapeNoMatchAtS/NMshapeMax * &
                                 (Langmuir_EFactor - cvmix_one)
      case default
         MixingCoefEnhancement = cvmix_one
      end select
      !   (3d) Diffusivity = OBL_depth * (turbulent scale) * G(sigma)
      OBL_Mdiff(kw) = OBL_depth * w_m(kw) * MshapeAtS * MixingCoefEnhancement
      OBL_Tdiff(kw) = OBL_depth * w_s(kw) * TshapeAtS * MixingCoefEnhancement
      OBL_Sdiff(kw) = OBL_depth * w_s(kw) * SshapeAtS * MixingCoefEnhancement
    end do

    ! (4) Compute the enhanced diffusivity
    !     (4a) Compute shape function at last cell center in OBL
    sigma_ktup = -zt(ktup)/OBL_depth
    MshapeAtS = cvmix_math_evaluate_cubic(Mshape, sigma_ktup)
    TshapeAtS = cvmix_math_evaluate_cubic(Tshape, sigma_ktup)
    SshapeAtS = cvmix_math_evaluate_cubic(Sshape, sigma_ktup)
    !     (4b) Compute turbulent scales at last cell center in OBL
    call cvmix_kpp_compute_turbulent_scales(sigma_ktup, OBL_depth, surf_buoy, &
                                       surf_fric, XIone, wm_ktup, ws_ktup, &
                                       CVmix_kpp_params_user)
    if (CVMix_KPP_Params_in%Langmuir_Mixing_Opt &
       .eq. LANGMUIR_MIXING_LWF16) then
      ! enhance the turbulent velocity scale
      wm_ktup = wm_ktup * Langmuir_EFactor
      ws_ktup = ws_ktup * Langmuir_EFactor
    end if
    !     (4c) Diffusivity = OBL_depth * (turbulent scale) * G(sigma)
    Mdiff_ktup = OBL_depth * wm_ktup * MshapeAtS
    Tdiff_ktup = OBL_depth * ws_ktup * TshapeAtS
    Sdiff_ktup = OBL_depth * ws_ktup * SshapeAtS

    if (CVmix_kpp_params_in%lenhanced_diff) then
      if ((ktup.eq.kwup).or.(ktup.eq.kwup-1)) then
        call cvmix_kpp_compute_enhanced_diff(Mdiff_ktup,                      &
                                             Tdiff_ktup,                      &
                                             Sdiff_ktup,                      &
                                             Mdiff_out(ktup+1),               &
                                             Tdiff_out(ktup+1),               &
                                             Sdiff_out(ktup+1),               &
                                             OBL_Mdiff(ktup+1),               &
                                             OBL_Tdiff(ktup+1),               &
                                             OBL_Sdiff(ktup+1),               &
                                             Tnonlocal(ktup+1),               &
                                             Snonlocal(ktup+1),               &
                                             delta, lkteqkw=(ktup.eq.kwup))
      else
        print*, "ERROR: ktup should be either kwup or kwup-1!"
        print*, "ktup = ", ktup, " and kwup = ", kwup
        stop 1
      end if
    else
      if ( kwup .eq. ktup ) then
        OBL_Mdiff(ktup+1) = old_Mdiff(ktup+1)
        OBL_Tdiff(ktup+1) = old_Tdiff(ktup+1)
        OBL_Sdiff(ktup+1) = old_Sdiff(ktup+1)
      end if
    end if

    ! (5) Combine interior and boundary coefficients
    Mdiff_out(2:ktup+1) = OBL_Mdiff(2:ktup+1)
    Tdiff_out(2:ktup+1) = OBL_Tdiff(2:ktup+1)
    Sdiff_out(2:ktup+1) = OBL_Sdiff(2:ktup+1)

  end if    ! lStokesMOST

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
      case ('minOBLdepth')
        CVmix_kpp_params_out%minOBLdepth = val
      case ('maxOBLdepth')
        CVmix_kpp_params_out%maxOBLdepth = val
      case ('minVtsqr')
        CVmix_kpp_params_out%minVtsqr = val
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
      case ('Cv')
        CVmix_kpp_params_out%Cv = val
      case ('CVt2')
        CVmix_kpp_params_out%CVt2 = val
      case ('nonlocal_coeff')
        CVmix_kpp_params_out%nonlocal_coeff = val
      case ('c_CT')
        CVmix_kpp_params_out%c_CT = val
      case ('c_ST')
        CVmix_kpp_params_out%c_ST = val
      case ('c_LT')
        CVmix_kpp_params_out%c_LT = val
      case ('p_LT')
         CVmix_kpp_params_out%p_LT = val
      case ('RWHGK_ENTR_COEF')
         CVmix_kpp_params_out%rwhgk_entr_coef = val
      case ('RWHGK_ENTR_EXP')
         CVmix_kpp_params_out%rwhgk_entr_exp = val
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
    end select

!EOC

  end subroutine cvmix_put_kpp_real

!BOP

! !IROUTINE: cvmix_put_kpp_int
! !INTERFACE:

  subroutine cvmix_put_kpp_int(varname, val, CVmix_kpp_params_user)

! !DESCRIPTION:
!  Write an integer value into a cvmix\_kpp\_params\_type variable.
!\\
!\\

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    integer,          intent(in) :: val

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
      case ('interp_type')
        CVmix_kpp_params_out%interp_type = val
      case ('interp_type2')
        CVmix_kpp_params_out%interp_type2 = val
      case ('MatchTechnique')
        CVmix_kpp_params_out%MatchTechnique = val
      case ('old_vals', 'handle_old_vals')
        CVmix_kpp_params_out%handle_old_vals = val
      case ('Langmuir_Mixing_Opt')
        CVmix_kpp_params_out%Langmuir_Mixing_opt = val
      case ('Langmuir_Entrainment_Opt')
        CVmix_kpp_params_out%Langmuir_Entrainment_opt = val
      case DEFAULT
        call cvmix_put_kpp(varname, real(val, cvmix_r8), CVmix_kpp_params_out)
    end select

!EOC

  end subroutine cvmix_put_kpp_int

!BOP

! !IROUTINE: cvmix_put_kpp_logical
! !INTERFACE:

  subroutine cvmix_put_kpp_logical(varname, val, CVmix_kpp_params_user)

! !DESCRIPTION:
!  Write a Boolean value into a cvmix\_kpp\_params\_type variable.
!\\
!\\

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    logical,          intent(in) :: val

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
      case ('lscalar_Cv')
        CVmix_kpp_params_out%lscalar_Cv = val
      case ('lEkman')
        CVmix_kpp_params_out%lEkman = val
      case ('lStokesMOST')
        CVmix_kpp_params_out%lStokesMOST = val
      case ('lMonOb')
        CVmix_kpp_params_out%lMonOb = val
      case ('lnoDGat1')
        CVmix_kpp_params_out%lnoDGat1 = val
      case ('lenhanced_diff')
        CVmix_kpp_params_out%lenhanced_diff = val
      case ('l_LMD_ws')
        CVmix_kpp_params_out%l_LMD_ws = val
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

! !INPUT PARAMETERS:
    character(len=*),                              intent(in) :: varname
    type(cvmix_kpp_params_type), optional, target, intent(in) ::              &
                                           CVmix_kpp_params_user

! !OUTPUT PARAMETERS:
    real(cvmix_r8) :: cvmix_get_kpp_real

!EOP
!BOC

    type(cvmix_kpp_params_type), pointer :: CVmix_kpp_params_get

    CVmix_kpp_params_get => CVmix_kpp_params_saved
    if (present(CVmix_kpp_params_user)) then
      CVmix_kpp_params_get => CVmix_kpp_params_user
    end if

    cvmix_get_kpp_real = cvmix_zero
    select case (trim(varname))
      case ('Ri_crit')
        cvmix_get_kpp_real = CVmix_kpp_params_get%Ri_crit
      case ('vonkarman')
        cvmix_get_kpp_real = CVmix_kpp_params_get%vonkarman
      case ('Cstar')
        cvmix_get_kpp_real = CVmix_kpp_params_get%Cstar
      case ('zeta_m')
        cvmix_get_kpp_real = CVmix_kpp_params_get%zeta_m
      case ('zeta_s')
        cvmix_get_kpp_real = CVmix_kpp_params_get%zeta_s
      case ('a_m')
        cvmix_get_kpp_real = CVmix_kpp_params_get%a_m
      case ('a_s')
        cvmix_get_kpp_real = CVmix_kpp_params_get%a_s
      case ('c_m')
        cvmix_get_kpp_real = CVmix_kpp_params_get%c_m
      case ('c_s')
        cvmix_get_kpp_real = CVmix_kpp_params_get%c_s
      case ('surf_layer_ext')
        cvmix_get_kpp_real = CVmix_kpp_params_get%surf_layer_ext
      case ('Cv')
        cvmix_get_kpp_real = CVmix_kpp_params_get%Cv
      case ('CVt2')
        cvmix_get_kpp_real = CVmix_kpp_params_get%CVt2
      case ('c_CT')
        cvmix_get_kpp_real = CVmix_kpp_params_get%c_CT
      case ('c_ST')
        cvmix_get_kpp_real = CVmix_kpp_params_get%c_ST
      case ('c_LT')
        cvmix_get_kpp_real = CVmix_kpp_params_get%c_LT
      case ('p_LT')
        cvmix_get_kpp_real = CVmix_kpp_params_get%p_LT
      case ('RWHGK_ENTR_COEF')
         cvmix_get_kpp_real = CVmix_kpp_params_get%RWHGK_ENTR_COEF
      case ('RWHGK_ENTR_EXP')
         cvmix_get_kpp_real = CVmix_kpp_params_get%RWHGK_ENTR_EXP
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
                                             surf_buoy, Coriolis,Xi, zBottom, &
                                             CVmix_kpp_params_user)

! !DESCRIPTION:
!  Computes the depth of the ocean boundary layer (\verb|OBL_depth|) for a given column.
!  Ri\_bulk(h) = Ricr; h < -zBottom, (stable+lMonOb) 0 < h < vonKaraman Lstar
!\\
!\\

! !INPUT PARAMETERS:
    real(cvmix_r8), dimension(:),                   intent(in) :: Ri_bulk
    real(cvmix_r8), dimension(:),           target, intent(in) :: zw_iface,   &
                                                                  zt_cntr
    real(cvmix_r8), dimension(:), optional,  intent(in) ::  Xi, surf_buoy
    real(cvmix_r8),               optional,         intent(in) :: surf_fric,  &
                                                                  Coriolis,   &
                                                                  zBottom
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
    integer             :: nlev, k, kRi
    logical             :: lstable

    type(cvmix_kpp_params_type), pointer :: CVmix_kpp_params_in

    CVmix_kpp_params_in => CVmix_kpp_params_saved
    if (present(CVmix_kpp_params_user)) then
      CVmix_kpp_params_in => CVmix_kpp_params_user
    end if

    ! Error checks
    ! (1) if using Ekman length, need to pass surf_fric and Coriolis
    if ((.not.(present(surf_fric).and.present(Coriolis))).and.                &
        CVmix_kpp_params_in%lEkman) then
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

    ! (3) zt_cntr must be length nlev and zw_iface must be nlev+1
    nlev = size(zt_cntr)
    if (size(zw_iface).ne.nlev+1) then
      print*, "ERROR: zt_cntr must have exactly one less element than zw_iface!"
      print*, "size(zt_cntr) = ", nlev, ", size(zw_iface) = ", size(zw_iface)
      stop 1
    end if

    ! (4) Ri_bulk needs to be either the size of zw_iface or zt_cntr
    if (size(Ri_bulk).eq.nlev) then
      depth => zt_cntr
    else if (size(Ri_bulk).eq.nlev+1) then
      depth => zw_iface
    else
      print*, "ERROR: Ri_bulk must have size nlev or nlev+1!"
      print*, "nlev = ", nlev, ", size(Ri_bulk) = ", size(Ri_bulk)
      stop 1
    end if

    ! if lEkman = .true., OBL_depth must be between the surface and the Ekman
    ! depth. Similarly, if lMonOb = .true., OBL_depth must be between the
    ! surface and the Monin-Obukhov depth

  if ( CVmix_kpp_params_in%lStokesMOST ) then
    ! OBL_depth must be at or above 1) zbottom, the effective ocean bottom,
    if ( present(zBottom) )  then
      OBL_limit = abs(zBottom)
    else
      OBL_limit = abs(zt_cntr(nlev))
    end if

    ! (1) Find k such that Ri_bulk at level k+1 > Ri_crit
    do k=0,size(Ri_bulk)-1
      kRi = k+1
      if (Ri_bulk(k+1).gt.CVmix_kpp_params_in%ri_crit) &
      exit
    end do

    if (k.eq.size(Ri_bulk)) then
      OBL_depth = OBL_limit
    elseif (k.eq.0) then
      OBL_depth = abs(zt_cntr(1))
    else

    ! (2) Interpolation
      if (k.eq.1) then
        call cvmix_math_poly_interp(coeffs, CVmix_kpp_params_in%interp_type,  &
                               depth(k:k+1), Ri_bulk(k:k+1))
      else
        call cvmix_math_poly_interp(coeffs, CVmix_kpp_params_in%interp_type,  &
                               depth(k:k+1), Ri_bulk(k:k+1), depth(k-1),      &
                               Ri_bulk(k-1))
      end if
      coeffs(1) = coeffs(1)-CVmix_kpp_params_in%ri_crit

      OBL_depth = -cvmix_math_cubic_root_find(coeffs, 0.5_cvmix_r8 *          &
                                                      (depth(k)+depth(k+1)))

    ! (3) OBL_depth needs to be at or below the center of the top level
      ! Note: OBL_depth can only be computed to be above this point if k=1,
      if (k.eq.1)  OBL_depth = max(OBL_depth, -zt_cntr(1))
    end if                 !  -zt_cntr(1)  <  OBL_depth < OBL_limit

      ! (5) the modified MoninObukhov limit (= vonk*Lstar) if stable and lMonOb=True
    if (CVmix_kpp_params_in%lMonOb ) then
      if ( present(Xi) .and. present(surf_buoy) ) then
        MoninObukhov = OBL_limit
        do k = 0, kRi-1
          if (surf_buoy(k+1) .gt. cvmix_zero) MoninObukhov =      &
                surf_fric**3 / (surf_buoy(k+1) * (cvmix_one-Xi(k+1)))
          if ( MoninObukhov .lt. abs(zt_cntr(k+1)) ) &
          exit
        end do
        if (k.eq.0) then
          OBL_limit  = abs(zt_cntr(1))
        elseif (k.lt.kRi) then
           OBL_limit = min( OBL_limit, abs(zw_iface(k)) )
        end if

      else
        print*, "ERROR: Stokes_XI and surf_buoy both must be present if lMonOb=true with Stokes_Most package"
        stop 1
      end if
    end if     ! lMonOb

      ! (4) OBL_depth must be at or above OBL_limit -zt_cntr(1)  <  OBL_depth < OBL_limit
    OBL_depth  = min(OBL_depth, OBL_limit )
    kOBL_depth = cvmix_kpp_compute_kOBL_depth(zw_iface, zt_cntr, OBL_depth)

  else ! not Stokes_MOST

    OBL_limit  = abs(zt_cntr(nlev))

    ! Since depth gets more negative as you go deeper, that translates into
    ! OBL_depth = max(abs(computed depth), abs(Ekman depth), abs(M-O depth))
    if (CVmix_kpp_params_in%lEkman) then
      ! Column is stable if surf_buoy > 0
      lstable = (surf_buoy(nlev).gt.cvmix_zero)

      if (Coriolis.ne.cvmix_zero .and. lstable) then
        Ekman = 0.7_cvmix_r8*surf_fric/abs(Coriolis)
      else
        ! Rather than divide by zero (or if column is unstable), set Ekman depth to ocean bottom
        Ekman = abs(zt_cntr(nlev))
      end if
      OBL_limit = min(OBL_limit, Ekman)
    end if

    if (CVmix_kpp_params_in%lMonOb) then
      ! Column is stable if surf_buoy > 0
      lstable = (surf_buoy(nlev).gt.cvmix_zero)

      if (lstable) then
        MoninObukhov = surf_fric**3/(surf_buoy(nlev)*CVmix_kpp_params_in%vonkarman)
      else
        MoninObukhov = abs(zt_cntr(nlev))
      end if
      OBL_limit = min(OBL_limit, MoninObukhov)
    end if

    ! Interpolation Step
    ! (1) Find k such that Ri_bulk at level k+1 > Ri_crit
    do k=0,size(Ri_bulk)-1
      if (Ri_bulk(k+1).gt.CVmix_kpp_params_in%ri_crit) &
        exit
    end do

    if (k.eq.size(Ri_bulk)) then
      OBL_depth = abs(OBL_limit)
    elseif (k.eq.0) then
      OBL_depth = abs(zt_cntr(1))
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

      OBL_depth = -cvmix_math_cubic_root_find(coeffs, 0.5_cvmix_r8 *          &
                                                      (depth(k)+depth(k+1)))

      ! OBL_depth needs to be at or below the center of the top level
      ! Note: OBL_depth can only be computed to be above this point if k=1,
      !       depth => zw_iface instead of zt_cntr, and the interpolation
      !       results in Ri_bulk = Ri_crit at a depth above the center of the
      !       top level.
      if (k.eq.1) then
        OBL_depth = max(OBL_depth, -zt_cntr(1))
      end if

      ! OBL_depth needs to be at or above OBL_limit
      ! Note: maybe there are times when we don't need to do the interpolation
      !       because we know OBL_depth will equal OBL_limit?
      OBL_depth = min(OBL_depth, OBL_limit)
    end if

    OBL_depth = max(OBL_depth, CVmix_kpp_params_in%minOBLdepth)
    if (CVmix_kpp_params_in%maxOBLdepth.gt.cvmix_zero)                        &
      OBL_depth = min(OBL_depth, CVmix_kpp_params_in%maxOBLdepth)
    kOBL_depth = cvmix_kpp_compute_kOBL_depth(zw_iface, zt_cntr, OBL_depth)

  end if    ! lStokesMOST

!EOC

  end subroutine cvmix_kpp_compute_OBL_depth_low

!BOP

! !IROUTINE: cvmix_kpp_compute_kOBL_depth
! !INTERFACE:

  function cvmix_kpp_compute_kOBL_depth(zw_iface, zt_cntr, OBL_depth)

! !DESCRIPTION:
!  Computes the index of the level and interface above \verb|OBL_depth|. The index is
!  stored as a real number, and the integer index can be solved for in the
!  following way:\\
!    \verb|kt| = index of cell center above \verb|OBL_depth| = \verb|nint(kOBL_depth)-1|
!    \verb|kw| = index of interface above \verb|OBL_depth| = \verb|floor(kOBL_depth)|
!\\
!\\

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

! !IROUTINE: cvmix_kpp_compute_enhanced_diff
! !INTERFACE:

  subroutine cvmix_kpp_compute_enhanced_diff(Mdiff_ktup, Tdiff_ktup,          &
                                             Sdiff_ktup, Mdiff, Tdiff, Sdiff, &
                                             OBL_Mdiff, OBL_Tdiff, OBL_Sdiff, &
                                             Tnonlocal, Snonlocal,            &
                                             delta, lkteqkw)

! !DESCRIPTION:
!  The enhanced mixing described in Appendix D of LMD94 changes the diffusivity
!  values at the interface between the cell center above \verb|OBL_depth| and the one
!  below it, based on a weighted average of how close to each center \verb|OBL_depth|
!  is. Note that we need to know whether \verb|OBL_depth| is above this interface or
!  below it - we do this by comparing the indexes of the cell center above
!  \verb|OBL_depth| (\verb|ktup|) and the cell interface above \verb|OBL_depth|(\verb|kwup|).
!\\
!\\

! !INPUT PARAMETERS:

    ! Diffusivity and viscosity at cell center above OBL_depth
    real(cvmix_r8), intent(in) :: Mdiff_ktup, Tdiff_ktup, Sdiff_ktup

    ! Weight to use in averaging (distance between OBL_depth and cell center
    ! above OBL_depth divided by distance between cell centers bracketing
    ! OBL_depth).
    real(cvmix_r8), intent(in) :: delta

    logical, intent(in) :: lkteqkw ! .true.  => interface ktup+1 is outside OBL
                                   !            (update diff and visc)
                                   ! .false. => interface ktup+1 is inside OBL
                                   !            (update OBL_diff and OBL_visc)

! !OUTPUT PARAMETERS:
    ! Will change either diff & visc or OBL_diff & OBL_visc, depending on value
    ! of lkteqkw
    real(cvmix_r8), intent(inout) :: Mdiff, Tdiff, Sdiff,                     &
                                     OBL_Mdiff, OBL_Tdiff, OBL_Sdiff,         &
                                     Tnonlocal, Snonlocal

!EOP
!BOC

    ! Local variables

    ! enh_diff and enh_visc are the enhanced diffusivity and viscosity values
    ! at the interface nearest OBL_depth
    real(cvmix_r8) :: enh_Mdiff, enh_Tdiff, enh_Sdiff
    ! Need to store original OBL_Tdiff and OBL_Sdiff for updating nonlocal
    real(cvmix_r8) :: old_Tdiff, old_Sdiff

    real(cvmix_r8) :: omd ! one minus delta

    omd = cvmix_one - delta
    old_Tdiff = OBL_Tdiff
    old_Sdiff = OBL_Sdiff

    if (lkteqkw) then
      ! => ktup = kwup
      ! Interface kw = ktup+1 is outside the OBL

      ! (a) compute enhanced diffs: get diffusivity values at kw = ktup+1
      !     from diff and visc rather than OBL_diff and OBL_visc
      enh_Mdiff = (omd**2)*Mdiff_ktup + (delta**2)*Mdiff
      enh_Tdiff = (omd**2)*Tdiff_ktup + (delta**2)*Tdiff
      enh_Sdiff = (omd**2)*Sdiff_ktup + (delta**2)*Sdiff

      ! (b) modify diffusivity values at kw = ktup+1 (again in diff and visc)
      Mdiff = omd*Mdiff + delta*enh_Mdiff
      Tdiff = omd*Tdiff + delta*enh_Tdiff
      Sdiff = omd*Sdiff + delta*enh_Sdiff

      ! (c) Update OBL_[MTS]diff
      OBL_Mdiff = Mdiff
      OBL_Tdiff = Tdiff
      OBL_Sdiff = Sdiff
    else
      ! => ktup = kwup - 1
      ! Interface kw = ktup+1 is in the OBL

      ! (a) compute enhanced diffs: get diffusivity values at kw = ktup+1
      !     from OBL_diff and OBL_visc rather than diff and visc
      enh_Mdiff = (omd**2)*Mdiff_ktup + (delta**2)*OBL_Mdiff
      enh_Tdiff = (omd**2)*Tdiff_ktup + (delta**2)*OBL_Tdiff
      enh_Sdiff = (omd**2)*Sdiff_ktup + (delta**2)*OBL_Sdiff

      ! (b) modify diffusivity values at kw = ktup+1 (again in OBL_diff and
      !     OBL_visc)
      OBL_Mdiff = omd*Mdiff + delta*enh_Mdiff
      OBL_Tdiff = omd*Tdiff + delta*enh_Tdiff
      OBL_Sdiff = omd*Sdiff + delta*enh_Sdiff

      ! (c) update nonlocal term
      if (old_Tdiff.ne.cvmix_zero) then
        Tnonlocal = Tnonlocal*OBL_Tdiff/old_Tdiff
      else
        Tnonlocal = cvmix_zero
      end if
      if (old_Sdiff.ne.cvmix_zero) then
        Snonlocal = Snonlocal*OBL_Sdiff/old_Sdiff
      else
        Snonlocal = cvmix_zero
      end if
    end if

! EOC

  end subroutine cvmix_kpp_compute_enhanced_diff

!BOP

! !IROUTINE: cvmix_kpp_compute_OBL_depth_wrap
! !INTERFACE:

  subroutine cvmix_kpp_compute_OBL_depth_wrap(CVmix_vars, CVmix_kpp_params_user)

! !DESCRIPTION:
!  Computes the depth of the ocean boundary layer (\verb|CVmix_vars%BoundaryLayerDepth|)
!  for a given column.
!\\
!\\

! !INPUT PARAMETERS:
    type(cvmix_kpp_params_type), optional, target, intent(in) ::                &
                                           CVmix_kpp_params_user

! !OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars

!EOP
!BOC

    ! Local variables
!   real(cvmix_r8) :: lcl_obl_depth, lcl_kobl_depth
    real(cvmix_r8), dimension(1) :: lcl_Xi, lcl_sfcBuoy
    lcl_Xi(1) = CVmix_vars%StokesMostXi
    lcl_sfcBuoy(1) = CVmix_vars%SurfaceBuoyancyForcing

    call cvmix_kpp_compute_OBL_depth(CVmix_vars%BulkRichardson_cntr,          &
                                     CVmix_vars%zw_iface,                     &
                                     CVmix_vars%BoundaryLayerDepth,           &
                                     CVmix_vars%kOBL_depth,                   &
                                     CVmix_vars%zt_cntr,                      &
                                     CVmix_vars%SurfaceFriction,              &
                                     lcl_sfcBuoy,                             &
                                     CVmix_vars%Coriolis,                     &
                                     lcl_Xi   ,                               &
                                     CVmix_vars%zBottomOceanNumerics,         &
                                     CVmix_kpp_params_user)

!EOC

  end subroutine cvmix_kpp_compute_OBL_depth_wrap

!BOP

! !IROUTINE: cvmix_kpp_compute_bulk_Richardson
! !INTERFACE:

  function cvmix_kpp_compute_bulk_Richardson(zt_cntr, delta_buoy_cntr,        &
                                             delta_Vsqr_cntr, Vt_sqr_cntr,    &
                                             ws_cntr, N_iface, Nsqr_iface,    &
                                             EFactor, LaSL, bfsfc, ustar,     &
                                             CVmix_kpp_params_user)

! !DESCRIPTION:
!  Computes the bulk Richardson number at cell centers. If \verb|Vt_sqr_cntr|
!  is not present, this routine will call \verb|compute_unresolved_shear|,
!  a routine that requires \verb|ws_cntr| and either \verb|N_iface| or
!  \verb|Nsqr_iface|.
!\\
!\\

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
    ! * bfsfc: surface buoyancy flux (units: m^2/s^3)
    real(cvmix_r8), dimension(size(zt_cntr)), intent(in), optional ::         &
                                              bfsfc, ws_cntr, Vt_sqr_cntr
    real(cvmix_r8), dimension(size(zt_cntr)+1), intent(in), optional ::       &
                                                    N_iface, Nsqr_iface
    ! * EFactor: Langmuir enhancement factor for entrainment (units: none)
    ! * LaSL: surface layer averaged Langmuir number (units: none)
    ! * ustar: friction velocity (units: m/s)
    real(cvmix_r8), intent(in), optional :: EFactor, LaSL, ustar
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
    real(cvmix_r8), dimension(size(zt_cntr)) :: unresolved_shear_cntr_sqr
    integer        :: kt
    real(cvmix_r8) :: scaling, num, denom
    type(cvmix_kpp_params_type),  pointer :: CVmix_kpp_params_in

    CVmix_kpp_params_in => CVmix_kpp_params_saved
    if (present(CVmix_kpp_params_user)) then
      CVmix_kpp_params_in => CVmix_kpp_params_user
    end if

    ! Make sure all arguments are same size
    if (any((/size(delta_buoy_cntr), size(delta_Vsqr_cntr)/).ne.              &
        size(zt_cntr))) then
      print*, "ERROR: delta_buoy, delta_vel_sqr, and zt_cntr must all be the",&
              "same size!"
      stop 1
    end if
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
      unresolved_shear_cntr_sqr = cvmix_kpp_compute_unresolved_shear(         &
                                    zt_cntr, ws_cntr, N_iface, Nsqr_iface,    &
                                    EFactor, LaSL, bfsfc, ustar,              &
                                    CVmix_kpp_params_user)
    end if

    ! scaling because we want (d-dr) = (d-0.5*eps*d) = (1-0.5*eps)*d
    scaling = cvmix_one - 0.5_cvmix_r8*CVmix_kpp_params_in%surf_layer_ext
    do kt=1,size(zt_cntr)
      ! Negative sign because we use positive-up for height
      num   = -scaling*zt_cntr(kt)*delta_buoy_cntr(kt)
      denom = delta_Vsqr_cntr(kt) + unresolved_shear_cntr_sqr(kt)
      if (denom.ne.cvmix_zero) then
        cvmix_kpp_compute_bulk_Richardson(kt) = num/denom
      else
        ! Need a better fudge factor?
        cvmix_kpp_compute_bulk_Richardson(kt) = num*1e10_cvmix_r8
      end if
    end do

!EOC

  end function cvmix_kpp_compute_bulk_Richardson

!BOP

! !IROUTINE: cvmix_kpp_compute_turbulent_scales_0d
! !INTERFACE:

  subroutine cvmix_kpp_compute_turbulent_scales_0d(sigma_coord, OBL_depth,    &
                                                   surf_buoy_force,           &
                                                   surf_fric_vel,             &
                                                   xi, w_m, w_s,              &
                                                   CVmix_kpp_params_user)

! !DESCRIPTION:
!  Computes the turbulent velocity scales for momentum (\verb|w_m|) and scalars
!  (\verb|w_s|) at a single $\sigma$ coordinate.
!\\
!\\

! !INPUT PARAMETERS:
    real(cvmix_r8), intent(in) :: sigma_coord
    real(cvmix_r8), intent(in) :: OBL_depth, surf_buoy_force, surf_fric_vel
    real(cvmix_r8), optional, intent(in) :: xi
    type(cvmix_kpp_params_type), intent(in), optional, target ::              &
                                           CVmix_kpp_params_user

! !OUTPUT PARAMETERS:
    real(cvmix_r8), optional, intent(inout) :: w_m
    real(cvmix_r8), optional, intent(inout) :: w_s

!EOP
!BOC

    ! Local variables
    real(cvmix_r8), dimension(1) :: sigma, lcl_wm, lcl_ws
    real(cvmix_r8)               :: lcl_XI
    logical :: compute_wm, compute_ws

    if( present( xi ) ) then
      lcl_XI = xi
    else
      lcl_XI = cvmix_zero           ! NO WAVES
    end if

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
                                              xi  = lcl_XI,                   &
                                              w_m = lcl_wm, w_s = lcl_ws,     &
                                  CVmix_kpp_params_user=CVmix_kpp_params_user)
    else
      if (compute_wm) &
        call cvmix_kpp_compute_turbulent_scales(sigma, OBL_depth,             &
                                                surf_buoy_force,surf_fric_vel,&
                                                xi  = lcl_XI,                 &
                                                w_m = lcl_wm,                 &
                                  CVmix_kpp_params_user=CVmix_kpp_params_user)
      if (compute_ws) &
        call cvmix_kpp_compute_turbulent_scales(sigma, OBL_depth,             &
                                                surf_buoy_force,surf_fric_vel,&
                                                xi  = lcl_XI,                 &
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

  subroutine cvmix_kpp_compute_turbulent_scales_1d_sigma(sigma_coord,         &
                                                         OBL_depth,           &
                                                         surf_buoy_force,     &
                                                         surf_fric_vel, xi,   &
                                                         w_m, w_s,            &
                                                         CVmix_kpp_params_user)

! !DESCRIPTION:
!  Computes the turbulent velocity scales for momentum (\verb|w_m|) and scalars
!  (\verb|w_s|) given a 1d array of $\sigma$ coordinates. Note that the
!  turbulent scales are a continuous function, so there is no restriction to
!  only evaluating this routine at interfaces or cell centers. Also, if
!  $\sigma >$ \verb|surf_layer_ext| (which is typically 0.1), \verb|w_m| and
!  \verb|w_s| will be evaluated at the latter value.
!\\
!\\

! !INPUT PARAMETERS:
    real(cvmix_r8), dimension(:), intent(in) :: sigma_coord
    real(cvmix_r8), intent(in) :: OBL_depth, surf_buoy_force, surf_fric_vel
    real(cvmix_r8), optional, intent(in)           :: xi
    type(cvmix_kpp_params_type), intent(in), optional, target ::              &
                                           CVmix_kpp_params_user

! !OUTPUT PARAMETERS:
    real(cvmix_r8), optional, dimension(size(sigma_coord)), intent(inout) ::  &
                                                                    w_m, w_s

!EOP
!BOC

    ! Local variables
    integer :: n_sigma, kw
    logical :: compute_wm, compute_ws, l_LMD_ws
    real(cvmix_r8), dimension(size(sigma_coord)) :: zeta, sigma_loc
    real(cvmix_r8) :: vonkar, surf_layer_ext
    real(cvmix_r8) :: chi_m, chi_s, L_StokesL
    type(cvmix_kpp_params_type), pointer :: CVmix_kpp_params_in

    n_sigma = size(sigma_coord)

    CVmix_kpp_params_in => CVmix_kpp_params_saved
    if (present(CVmix_kpp_params_user)) then
      CVmix_kpp_params_in => CVmix_kpp_params_user
    end if

    compute_wm = present(w_m)
    compute_ws = present(w_s)

    l_LMD_ws       = CVmix_kpp_params_in%l_LMD_ws
    vonkar         = CVmix_kpp_params_in%vonkarman
    surf_layer_ext = CVmix_kpp_params_in%surf_layer_ext

  if ( CVmix_kpp_params_in%lStokesMOST ) then
    if (present(xi)) then
      L_StokesL = cvmix_one - xi

      if (surf_fric_vel.gt.cvmix_zero) then
        sigma_loc(:) = min(cvmix_one , sigma_coord(:))
        zeta(:) = sigma_loc(:) * OBL_depth * surf_buoy_force * vonkar / &
              (surf_fric_vel**3)
        if (compute_wm) then
          chi_m = compute_Stokes_chi( xi , lchi_m=.true. )
          do kw=1,n_sigma
            w_m(kw) = compute_phi_inv(zeta(kw),CVmix_kpp_params_in, L_StokesL, lphi_m=.true.)*&
                      vonkar*surf_fric_vel / chi_m
          end do
        end if
        if (compute_ws) then
          chi_s = compute_Stokes_chi( xi , lchi_s=.true. )
          do kw=1,n_sigma
            w_s(kw) = compute_phi_inv(zeta(kw),CVmix_kpp_params_in, L_StokesL, lphi_s=.true.)*&
                      vonkar*surf_fric_vel / chi_s
          end do
        end if

      else ! surf_fric_vel = 0
        if (compute_wm) then
          if (surf_buoy_force.ge.cvmix_zero) then  ! STABLE
            w_m = cvmix_zero
          else                                 ! convective limit
            chi_m = compute_Stokes_chi( xi , lchi_m=.true. )
            L_StokesL = cvmix_one - xi
            do kw=1,n_sigma
              w_m(kw) = -surf_buoy_force * real(14,cvmix_r8) * sigma_coord(kw) * &
                         OBL_depth * vonkar * L_StokesL
              w_m(kw) = vonkar*(w_m(kw)**(cvmix_one/real(3,cvmix_r8))) / chi_m
            end do
          end if
        end if   ! compute_wm

        if (compute_ws) then
          if (surf_buoy_force.ge.cvmix_zero) then  ! STABLE
            w_s = cvmix_zero
          else
            chi_s = compute_Stokes_chi( xi , lchi_s=.true. )
            L_StokesL = cvmix_one - xi
            do kw=1,n_sigma                    ! convective limit
              w_s(kw) = -surf_buoy_force * real(25,cvmix_r8) * sigma_coord(kw) * &
                        OBL_depth * vonkar * L_StokesL
              w_s(kw) = vonkar*(w_s(kw)**(cvmix_one/real(3,cvmix_r8))) / chi_s
            end do
          end if  ! surf_buoy_force >= 0
        end if    ! compute_ws

      end if ! surf_fric_vel != 0

    else
      print*, "ERROR: Similarity xi must be present in 1d_sigma to use Stokes_MOST package!"
      stop 1
    end if    ! lStokesMOST and xi not present

  else    ! not lStokesMOST

    if (surf_fric_vel.ne.cvmix_zero) then
      if ((surf_buoy_force.ge.cvmix_zero) .and. l_LMD_ws) then
        sigma_loc(:) = sigma_coord(:)
      else
        sigma_loc(:) = min(surf_layer_ext, sigma_coord(:))
      end if
      ! compute scales at sigma if sigma < surf_layer_ext, otherwise compute
      ! at surf_layer_ext
      zeta(:) = sigma_loc(:) * OBL_depth * surf_buoy_force * vonkar /         &
                (surf_fric_vel**3)

      if (compute_wm) then
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
    else ! surf_fric_vel = 0
      if (compute_wm) then
        if (surf_buoy_force.ge.cvmix_zero) then
          ! Stable regime with surf_fric_vel = 0 => w_m = 0
          w_m = cvmix_zero
        else
          ! Unstable forcing, Eqs. (13) and (B1c) reduce to following
          do kw=1,n_sigma
            ! Compute (u*/phi_m)^3 [this is where the zeros in numerator and
            !                       denominator cancel when u* = 0]
            w_m(kw) = -CVmix_kpp_params_in%c_m *                              &
                      min(surf_layer_ext, sigma_coord(kw)) * OBL_depth *      &
                      vonkar * surf_buoy_force
            ! w_m = vonkar * u* / phi_m
            !     = vonkar * ((u*/phi_m)^3)^1/3
            w_m(kw) = vonkar*(w_m(kw)**(cvmix_one/real(3,cvmix_r8)))
          end do
        end if ! surf_buoy_force >= 0
      end if ! compute_wm

      if (compute_ws) then
        if (surf_buoy_force.ge.cvmix_zero) then
          ! Stable regime with surf_fric_vel = 0 => w_s = 0
          w_s = cvmix_zero
        else
          ! Unstable forcing, Eqs. (13) and (B1e) reduce to following
          do kw=1,n_sigma
            ! Compute (u*/phi_s)^3 [this is where the zeros in numerator and
            !                       denominator cancel when u* = 0]
            w_s(kw) = -CVmix_kpp_params_in%c_s *                              &
                      min(surf_layer_ext, sigma_coord(kw)) * OBL_depth *      &
                      vonkar * surf_buoy_force
            ! w_s = vonkar * u* / phi_s
            !     = vonkar * ((u*/phi_s)^3)^1/3
            w_s(kw) = vonkar*(w_s(kw)**(cvmix_one/real(3,cvmix_r8)))
          end do
        end if ! surf_buoy_force >= 0
      end if ! compute_ws
    end if ! surf_fric_vel != 0
  end if    ! lStokesMOST
!EOC

  end subroutine cvmix_kpp_compute_turbulent_scales_1d_sigma

!BOP

! !IROUTINE: cvmix_kpp_compute_turbulent_scales_1d_OBL
! !INTERFACE:

  subroutine cvmix_kpp_compute_turbulent_scales_1d_OBL(sigma_coord,           &
                                                       OBL_depth,             &
                                                       surf_buoy_force,       &
                                                       surf_fric_vel,         &
                                                       xi, w_m, w_s,          &
                                                       CVmix_kpp_params_user)

! !DESCRIPTION:
!  Computes the turbulent velocity scales for momentum (\verb|w_m|) and scalars
!  (\verb|w_s|) given a single $\sigma$ coordinate and an array of boundary
!  layer depths. Note that the turbulent scales are a continuous function, so
!  they are evaluated at sigma\_coord * OBL\_depth(z) using surf\_buoy\_force(z)
!\\
!\\

! !INPUT PARAMETERS:
    real(cvmix_r8), intent(in) :: sigma_coord
    real(cvmix_r8), intent(in) :: surf_fric_vel
    real(cvmix_r8), dimension(:), intent(in) ::  surf_buoy_force, OBL_depth
    real(cvmix_r8), dimension(:), intent(in), optional :: xi
    type(cvmix_kpp_params_type), intent(in), optional, target ::              &
                                           CVmix_kpp_params_user

! !OUTPUT PARAMETERS:
    real(cvmix_r8), optional, dimension(size(surf_buoy_force)), intent(inout) &
                                                                :: w_m, w_s

!EOP
!BOC

    ! Local variables
    integer :: n_sigma, kw
    logical :: compute_wm, compute_ws, l_LMD_ws
    real(cvmix_r8), dimension(size(surf_buoy_force)) :: zeta, sigma_loc
    real(cvmix_r8) :: vonkar, surf_layer_ext
    real(cvmix_r8)  :: chi_m,chi_s,L_StokesL
    type(cvmix_kpp_params_type), pointer :: CVmix_kpp_params_in

    n_sigma = size(surf_buoy_force)

    CVmix_kpp_params_in => CVmix_kpp_params_saved
    if (present(CVmix_kpp_params_user)) then
      CVmix_kpp_params_in => CVmix_kpp_params_user
    end if

    compute_wm = present(w_m)
    compute_ws = present(w_s)

    l_LMD_ws       = CVmix_kpp_params_in%l_LMD_ws
    vonkar         = CVmix_kpp_params_in%vonkarman
    surf_layer_ext = CVmix_kpp_params_in%surf_layer_ext

  if ( CVmix_kpp_params_in%lStokesMOST ) then
    if (present(xi)) then

      if (surf_fric_vel.ne.cvmix_zero) then
        zeta(:) = sigma_coord*OBL_depth(:)*surf_buoy_force(:)*vonkar / &
                     (surf_fric_vel**3)
        if (compute_wm) then
          do kw = 1,n_sigma
            chi_m = compute_Stokes_chi( xi(kw) , lchi_m=.true. )
            L_StokesL = cvmix_one - xi(kw)             !wgl - xi
            w_m(kw)=compute_phi_inv(zeta(kw),CVmix_kpp_params_in,L_StokesL,lphi_m=.true.) * &
                      vonkar*surf_fric_vel / chi_m
          end do
        end if

        if (compute_ws) then
          do kw = 1,n_sigma
            chi_s = compute_Stokes_chi( xi(kw) , lchi_s=.true. )
            L_StokesL = cvmix_one - xi(kw)             !wgl - xi
            w_s(kw)=compute_phi_inv(zeta(kw),CVmix_kpp_params_in,L_StokesL,lphi_s=.true.) * &
                      vonkar*surf_fric_vel  / chi_s
          end do
        end if

      else ! surf_fric_vel = 0
        if (compute_wm) then
          do kw = 1,n_sigma
            if (surf_buoy_force(kw).ge.cvmix_zero) then  ! STABLE
              w_m(kw) = cvmix_zero
            else                            ! convective limit
              L_StokesL = cvmix_one - xi(kw)
              w_m(kw) = -surf_buoy_force(kw) * real(14,cvmix_r8) * sigma_coord * &
                        OBL_depth(kw) * vonkar * L_StokesL
              w_m(kw) =vonkar*(w_m(kw)**(cvmix_one/real(3,cvmix_r8))) / compute_Stokes_chi(xi(kw),lchi_m=.true.)
            end if
          end do
        end if
        if (compute_ws) then
          do kw = 1,n_sigma
            if (surf_buoy_force(kw).ge.cvmix_zero) then  ! STABLE
              w_s(kw) = cvmix_zero
            else
              L_StokesL = cvmix_one - xi(kw)
              w_s(kw) = -surf_buoy_force(kw) * real(25,cvmix_r8) * sigma_coord * &
                        OBL_depth(kw) * vonkar * L_StokesL
              w_s(kw) =vonkar*(w_s(kw)**(cvmix_one/real(3,cvmix_r8))) / compute_Stokes_chi(xi(kw),lchi_s=.true.)
            end if  ! surf_buoy_force >= 0
          end do
        end if    ! compute_ws
      end if      ! surf_fric_vel = 0
    else
      print*, "ERROR: Similarity xi must be present in 1d_OBL to use Stokes_MOST package!"
      stop 1
    end if    ! lStokesMOST but xi not present

  else    ! not lStokesMOST

    if (surf_fric_vel.ne.cvmix_zero) then
      sigma_loc = min(surf_layer_ext, sigma_coord)
      if (l_LMD_ws) then
        where (surf_buoy_force.ge.cvmix_zero)
          sigma_loc = sigma_coord
        end where
      end if
      zeta(:) = sigma_loc(:) * OBL_depth(:) * surf_buoy_force(:) * vonkar /   &
                (surf_fric_vel**3)

      if (compute_wm) then
        w_m(1) = compute_phi_inv(zeta(1), CVmix_kpp_params_in, lphi_m=.true.)*&
                 vonkar*surf_fric_vel
        do kw=2,n_sigma
          if (zeta(kw).eq.zeta(kw-1)) then
            w_m(kw) = w_m(kw-1)
          else
            w_m(kw) = compute_phi_inv(zeta(kw), CVmix_kpp_params_in, lphi_m=.true.)*&
                      vonkar*surf_fric_vel
          end if
        end do
      end if

      if (compute_ws) then
        w_s(1) = compute_phi_inv(zeta(1), CVmix_kpp_params_in, lphi_s=.true.)*&
                 vonkar*surf_fric_vel
        do kw=2,n_sigma
          if (zeta(kw).eq.zeta(kw-1)) then
            w_s(kw) = w_s(kw-1)
          else
            w_s(kw) = compute_phi_inv(zeta(kw), CVmix_kpp_params_in, lphi_s=.true.)*&
                      vonkar*surf_fric_vel
          end if
        end do
      end if

    else ! surf_fric_vel = 0
      if (compute_wm) then
        ! Unstable forcing, Eqs. (13) and (B1c) reduce to following
        do kw=1,n_sigma
          if(surf_buoy_force(kw) .ge. cvmix_zero) then
            w_m(kw) = cvmix_zero
          else
            ! Compute (u*/phi_m)^3 [this is where the zeros in numerator and
            !                       denominator cancel when u* = 0]
            w_m(kw) = -CVmix_kpp_params_in%c_m *                              &
                      min(surf_layer_ext, sigma_coord) * OBL_depth(kw) *      &
                      vonkar * surf_buoy_force(kw)
            ! w_m = vonkar * u* / phi_m
            !     = vonkar * ((u*/phi_m)^3)^1/3
            w_m(kw) = vonkar*(w_m(kw)**(cvmix_one/real(3,cvmix_r8)))
        end if
        end do
      end if ! compute_wm

      if (compute_ws) then
          ! Unstable forcing, Eqs. (13) and (B1e) reduce to following
        do kw=1,n_sigma
          if (surf_buoy_force(kw) .ge. cvmix_zero) then
            ! Stable regime with surf_fric_vel = 0 => w_s = 0
            w_s(kw) = cvmix_zero
          else
            ! Unstable forcing, Eqs. (13) and (B1e) reduce to following
            ! Compute (u*/phi_s)^3 [this is where the zeros in numerator and
            !                       denominator cancel when u* = 0]
            w_s(kw) = -CVmix_kpp_params_in%c_s *                              &
                      min(surf_layer_ext, sigma_coord) * OBL_depth(kw) *      &
                      vonkar * surf_buoy_force(kw)
            ! w_s = vonkar * u* / phi_s
            !     = vonkar * ((u*/phi_s)^3)^1/3
            w_s(kw) = vonkar*(w_s(kw)**(cvmix_one/real(3,cvmix_r8)))
          end if ! surf_buoy_force >= 0
        end do
      end if ! compute_ws
    end if ! surf_fric_vel != 0
  end if    ! lStokesMOST

!EOC

  end subroutine cvmix_kpp_compute_turbulent_scales_1d_OBL

!BOP

! !IROUTINE: cvmix_kpp_compute_unresolved_shear
! !INTERFACE:

  function cvmix_kpp_compute_unresolved_shear(zt_cntr, ws_cntr, N_iface,      &
                                            Nsqr_iface, EFactor,              &
                                            LaSL, bfsfc, ustar,               &
                                            CVmix_kpp_params_user)

! !DESCRIPTION:
!  Computes the square of the unresolved shear ($V_t^2$ in Eq. (23) of LMD94)
!  at cell centers. Note that you must provide either the buoyancy frequency
!  or its square at cell interfaces, this routine by default will use the
!  lower cell interface value as the cell center, but you can instead take
!  an average of the top and bottom interface values by setting
!  lavg\_N\_or\_Nsqr = .true. in cvmix\_kpp\_init(). If you pass in Nsqr then
!  negative values are assumed to be zero (default POP behavior).
!\\
!\\

! !INPUT PARAMETERS:
    ! zt_cntr: height at center of cell (units: m)
    ! ws_cntr: w_s (turbulent scale factor) at center of cell (units: m/s)
    real(cvmix_r8), dimension(:), intent(in) :: zt_cntr,  ws_cntr
    ! N_iface: buoyancy frequency at cell interfaces (units: 1/s)
    ! Nsqr_iface: squared buoyancy frequency at cell interfaces (units: 1/s^2)
    ! note that you must provide exactly one of these two inputs!
    real(cvmix_r8), dimension(size(zt_cntr)+1), intent(in), optional ::       &
                                                    N_iface, Nsqr_iface
    ! bfsfc: surface buoyancy flux above cell centers (units: m^2/s^3)
    real(cvmix_r8), dimension(size(zt_cntr)), intent(in), optional :: bfsfc
    ! EFactor: Langmuir enhancement factor (units: none)
    ! LaSL: surface layer averaged Langmuir number (units: none)
    ! ustar: friction velocity (units: m/s)
    real(cvmix_r8), intent(in), optional :: EFactor, LaSL, ustar
    type(cvmix_kpp_params_type),  intent(in), optional, target ::             &
                                           CVmix_kpp_params_user

! !OUTPUT PARAMETERS:
    real(cvmix_r8), dimension(size(zt_cntr)) ::                               &
                             cvmix_kpp_compute_unresolved_shear

!EOP
!BOC

    ! Local variables
    integer :: kt, nlev
    real(cvmix_r8) :: Cv, Vtc, beta, xbeta, ybeta   ! pure convection entrainment rule
    logical :: lwstar          ! use wstar rather than w_s
    real(cvmix_r8) :: wstar    ! convective velocity scale
    real(cvmix_r8) :: ws_wstar ! ratio in limit of pure convection
    ! N_cntr: buoyancy frequency at cell centers, derived from either N_iface
    !        or Nsqr_iface (units: 1/s)
    real(cvmix_r8), dimension(size(zt_cntr)) :: N_cntr
    ! c_CT, c_ST, c_LT, p_LT: parameters of Langmuir-enhanced entrainment
    !                         in Li and Fox-Kemper, 2017, JPO
    real(cvmix_r8) :: c_CT, c_ST, c_LT, p_LT
    ! RWHGK_ENTR_COEF, RWHGK_ENTR_EXP: parameters of Langmuir-enhanced
    !                         entrainment in Reichl et al., 2016, JPO
    real(cvmix_r8) :: RWHGK_ENTR_COEF, RWHGK_ENTR_EXP
    ! Vt2_Enhancement: enhancement factor for unresolved shear
    real(cvmix_r8) :: Vt2_Enhancement
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

    CVmix_kpp_params_in => CVmix_kpp_params_saved
    if (present(CVmix_kpp_params_user)) then
      CVmix_kpp_params_in => CVmix_kpp_params_user
    end if

    if (present(N_iface)) then
      if (size(N_iface).ne.(nlev+1)) then
        print*, "ERROR: N_iface must have one more element than zt_cntr"
        stop 1
      end if
      do kt=1,nlev
        N_cntr(kt) = N_iface(kt+1)
      end do
    else
      if (present(Nsqr_iface)) then
        if (size(Nsqr_iface).ne.(nlev+1)) then
          print*, "ERROR: Nsqr_iface must have one more element than zt_cntr"
          stop 1
        end if
        do kt=1,nlev
          N_cntr(kt)=sqrt(max(Nsqr_iface(kt+1),cvmix_zero))
        end do
      else
        print*, "ERROR: you must provide N_iface OR Nsqr_iface"
        stop 1
      end if
    end if

  if ( CVmix_kpp_params_in%lStokesMOST ) then
    lwstar   = .false.   !  .true.
    ws_wstar = CVmix_kpp_params_in%vonkarman * cvmix_one * real(25,cvmix_r8) ! *  &
    ws_wstar = CVmix_kpp_params_in%vonkarman * ws_wstar**(cvmix_one/real(3,cvmix_r8))

    do kt=1,nlev
      xbeta = MAX( cvmix_zero , MIN( (300._cvmix_r8-100._cvmix_r8) , (-zt_cntr(kt) - 100._cvmix_r8) ) )
      xbeta = xbeta / (300._cvmix_r8-100._cvmix_r8)
      ybeta = (cvmix_one - xbeta*xbeta)**2
      beta  = ybeta * 0.2_cvmix_r8
      Vtc = sqrt( beta *3.8409_cvmix_r8 /ws_wstar) /CVmix_kpp_params_in%Ri_crit
      if (lwstar ) then
        wstar = (MAX(0.0 , zt_cntr(kt) * bfsfc(kt) ))**(cvmix_one/real(3,cvmix_r8))
        cvmix_kpp_compute_unresolved_shear(kt) = &
             -zt_cntr(kt) * N_cntr(kt)  * CVmix_kpp_params_in%CVt2 * Vtc * wstar
      else
        cvmix_kpp_compute_unresolved_shear(kt) = &
             -zt_cntr(kt) * N_cntr(kt)  * CVmix_kpp_params_in%CVt2 * Vtc * ws_cntr(kt) / ws_wstar
      end if

      if (cvmix_kpp_compute_unresolved_shear(kt).lt.                       &
              CVmix_kpp_params_in%minVtsqr) then
            cvmix_kpp_compute_unresolved_shear(kt) = CVmix_kpp_params_in%minVtsqr
      end if
    enddo
  else    ! not lStokesMOST

    ! options for Langmuir enhanced entrainment
    select case (CVmix_kpp_params_in%Langmuir_Entrainment_Opt)

      case (LANGMUIR_ENTRAINMENT_LWF16)
        if (.not.(present(EFactor) )) then
           print*, "ERROR: you must pass in EFactor if ",&
                "Langmuir_entrainment_str .eq. 'LWF16'!"
           stop 1
        end if
        Vt2_Enhancement = EFactor

        ! From LMD 94, Vtc = sqrt(-beta_T/(c_s*eps))/kappa^2
        Vtc = sqrt(0.2_cvmix_r8/(cvmix_get_kpp_real('c_s', CVmix_kpp_params_in) * &
              cvmix_get_kpp_real('surf_layer_ext', CVmix_kpp_params_in))) / &
              (cvmix_get_kpp_real('vonkarman', CVmix_kpp_params_in)**2)

        do kt=1,nlev
          if (CVmix_kpp_params_in%lscalar_Cv) then
            Cv = cvmix_get_kpp_real('Cv', CVmix_kpp_params_in)
          else
            ! Cv computation comes from Danabasoglu et al., 2006
            if (N_cntr(kt).lt.0.002_cvmix_r8) then
              Cv = 2.1_cvmix_r8-real(200,cvmix_r8)*N_cntr(kt)
            else
              Cv = 1.7_cvmix_r8
            end if
          end if

          cvmix_kpp_compute_unresolved_shear(kt) = -Cv*Vtc*zt_cntr(kt)*          &
                                N_cntr(kt)*ws_cntr(kt)/                          &
                                CVmix_kpp_params_in%Ri_crit * Vt2_Enhancement
          if (cvmix_kpp_compute_unresolved_shear(kt).lt.                         &
              CVmix_kpp_params_in%minVtsqr) then
            cvmix_kpp_compute_unresolved_shear(kt) = CVmix_kpp_params_in%minVtsqr
          end if
        end do

      case (LANGMUIR_ENTRAINMENT_LF17)

        if (.not.(present(LaSL) .and. present(bfsfc) .and. present(ustar))) then
          print*, "ERROR: you must pass in LaSL, bfsfc and ustar if ",&
                "Langmuir_entrainment_str == 'LF17'!"
          stop 1
        end if
        ! only apply Langmuir enhanced entrainment under unstable condition
        if (bfsfc(1)<cvmix_zero) then
          ! (26) of Li and Fox-Kemper, 2017, JPO
          c_CT =  cvmix_get_kpp_real('c_CT', CVmix_kpp_params_in)
          c_ST =  cvmix_get_kpp_real('c_ST', CVmix_kpp_params_in)
          c_LT =  cvmix_get_kpp_real('c_LT', CVmix_kpp_params_in)
          p_LT =  cvmix_get_kpp_real('p_LT', CVmix_kpp_params_in)
          do kt=1,nlev
            if (CVmix_kpp_params_in%lscalar_Cv) then
              Cv = cvmix_get_kpp_real('Cv', CVmix_kpp_params_in)
            else
            ! Cv computation comes from Danabasoglu et al., 2006
              if (N_cntr(kt).lt.0.002_cvmix_r8) then
                Cv = 2.1_cvmix_r8-real(200,cvmix_r8)*N_cntr(kt)
              else
                Cv = 1.7_cvmix_r8
              end if
            end if
            Vtc = sqrt((c_CT*bfsfc(kt)*zt_cntr(kt) + c_ST*ustar**3 +             &
                     c_LT*ustar**3*LaSL**(-1.*p_LT))/ws_cntr(kt))
            cvmix_kpp_compute_unresolved_shear(kt) = -Cv*Vtc*zt_cntr(kt)*        &
                                   N_cntr(kt)/CVmix_kpp_params_in%Ri_crit
            if (cvmix_kpp_compute_unresolved_shear(kt).lt.                       &
                CVmix_kpp_params_in%minVtsqr) then
              cvmix_kpp_compute_unresolved_shear(kt) = CVmix_kpp_params_in%minVtsqr
            end if
          end do
        else
          ! From LMD 94, Vtc = sqrt(-beta_T/(c_s*eps))/kappa^2
          Vtc = sqrt(0.2_cvmix_r8/(cvmix_get_kpp_real('c_s', CVmix_kpp_params_in) * &
                cvmix_get_kpp_real('surf_layer_ext', CVmix_kpp_params_in))) / &
                (cvmix_get_kpp_real('vonkarman', CVmix_kpp_params_in)**2)

          do kt=1,nlev
            if (CVmix_kpp_params_in%lscalar_Cv) then
              Cv = cvmix_get_kpp_real('Cv', CVmix_kpp_params_in)
            else
              ! Cv computation comes from Danabasoglu et al., 2006
              if (N_cntr(kt).lt.0.002_cvmix_r8) then
                Cv = 2.1_cvmix_r8-real(200,cvmix_r8)*N_cntr(kt)
              else
                Cv = 1.7_cvmix_r8
              end if
            end if

            cvmix_kpp_compute_unresolved_shear(kt) = -Cv*Vtc*zt_cntr(kt) *       &
                              N_cntr(kt)*ws_cntr(kt)/CVmix_kpp_params_in%Ri_crit
            if (cvmix_kpp_compute_unresolved_shear(kt).lt.                       &
                CVmix_kpp_params_in%minVtsqr) then
              cvmix_kpp_compute_unresolved_shear(kt) = CVmix_kpp_params_in%minVtsqr
            end if
          end do
        end if

      case (LANGMUIR_ENTRAINMENT_RWHGK16)

        if (.not.(present(LaSL) )) then
           print*, "ERROR: you must pass in LaSL if ",&
                "Langmuir_entrainment_str == 'RWHGK16'!"
           stop 1
        end if
        RWHGK_ENTR_COEF =  cvmix_get_kpp_real('RWHGK_ENTR_COEF', &
             CVmix_kpp_params_in)
        RWHGK_ENTR_EXP =  cvmix_get_kpp_real('RWHGK_ENTR_EXP', &
             CVmix_kpp_params_in)
        Vt2_Enhancement = cvmix_one + RWHGK_ENTR_COEF * LASL**RWHGK_ENTR_EXP

        ! From LMD 94, Vtc = sqrt(-beta_T/(c_s*eps))/kappa^2
        Vtc = sqrt(0.2_cvmix_r8/(cvmix_get_kpp_real('c_s', CVmix_kpp_params_in) * &
              cvmix_get_kpp_real('surf_layer_ext', CVmix_kpp_params_in))) / &
              (cvmix_get_kpp_real('vonkarman', CVmix_kpp_params_in)**2)

        do kt=1,nlev
          if (CVmix_kpp_params_in%lscalar_Cv) then
            Cv = cvmix_get_kpp_real('Cv', CVmix_kpp_params_in)
          else
            ! Cv computation comes from Danabasoglu et al., 2006
            if (N_cntr(kt).lt.0.002_cvmix_r8) then
              Cv = 2.1_cvmix_r8-real(200,cvmix_r8)*N_cntr(kt)
            else
              Cv = 1.7_cvmix_r8
            end if
          end if

          cvmix_kpp_compute_unresolved_shear(kt) = -Cv*Vtc*zt_cntr(kt)*          &
                                N_cntr(kt)*ws_cntr(kt)/                          &
                                CVmix_kpp_params_in%Ri_crit * Vt2_Enhancement
          if (cvmix_kpp_compute_unresolved_shear(kt).lt.                         &
              CVmix_kpp_params_in%minVtsqr) then
            cvmix_kpp_compute_unresolved_shear(kt) = CVmix_kpp_params_in%minVtsqr
          end if
        end do

      case DEFAULT

        ! From LMD 94, Vtc = sqrt(-beta_T/(c_s*eps))/kappa^2
        Vtc = sqrt(0.2_cvmix_r8/(cvmix_get_kpp_real('c_s', CVmix_kpp_params_in) * &
              cvmix_get_kpp_real('surf_layer_ext', CVmix_kpp_params_in))) / &
              (cvmix_get_kpp_real('vonkarman', CVmix_kpp_params_in)**2)

        do kt=1,nlev
          if (CVmix_kpp_params_in%lscalar_Cv) then
            Cv = cvmix_get_kpp_real('Cv', CVmix_kpp_params_in)
          else
            ! Cv computation comes from Danabasoglu et al., 2006
            if (N_cntr(kt).lt.0.002_cvmix_r8) then
              Cv = 2.1_cvmix_r8-real(200,cvmix_r8)*N_cntr(kt)
            else
              Cv = 1.7_cvmix_r8
            end if
          end if

          cvmix_kpp_compute_unresolved_shear(kt) = -Cv*Vtc*zt_cntr(kt) *       &
                            N_cntr(kt)*ws_cntr(kt)/CVmix_kpp_params_in%Ri_crit
          if (cvmix_kpp_compute_unresolved_shear(kt).lt.                       &
              CVmix_kpp_params_in%minVtsqr) then
            cvmix_kpp_compute_unresolved_shear(kt) = CVmix_kpp_params_in%minVtsqr
          end if
        end do

    end select
  end if    ! lStokesMOST

!EOC

  end function cvmix_kpp_compute_unresolved_shear

!BOP

! !IROUTINE: compute_phi_inv
! !INTERFACE:

  function compute_phi_inv(zeta, CVmix_kpp_params_in, L_Lstokes, lphi_m, lphi_s)

! !DESCRIPTION:
!  Computes $\frac{1}{\phi_m}$ or $\frac{1}{\phi_s}$
!\\
!\\

! !INPUT PARAMETERS:
    real(cvmix_r8),              intent(in) :: zeta
    type(cvmix_kpp_params_type), intent(in) :: CVmix_kpp_params_in
    real(cvmix_r8), optional,    intent(in) :: L_Lstokes
    logical, optional,           intent(in) :: lphi_m, lphi_s

! !OUTPUT PARAMETERS:
    real(cvmix_r8) :: compute_phi_inv

!EOP
!BOC

    logical :: lm, ls
    real(cvmix_r8) :: zetastar

    ! If not specifying lphi_m or lphi_s, routine will error out, but
    ! initializing result to 0 removes warning about possibly returning an
    ! un-initialized value
    compute_phi_inv = cvmix_zero

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

    zetastar =  zeta
    if ( present(L_Lstokes) ) zetastar = zeta * L_Lstokes

    if ( CVmix_kpp_params_in%lStokesMOST ) then
      if (lm) then
        if (zeta.ge.cvmix_zero) then   ! STABLE
          compute_phi_inv = cvmix_one/(cvmix_one + real(14,cvmix_r8)*zetastar)
        else
          compute_phi_inv = &
          (cvmix_one - real(14,cvmix_r8)*zetastar)**(cvmix_one/real(3,cvmix_r8))
        end if
      end if
      if (ls) then
        if (zeta.ge.cvmix_zero) then   ! STABLE
          compute_phi_inv = cvmix_one/(cvmix_one + real( 5,cvmix_r8)*zetastar)
        else
          compute_phi_inv = &
          (cvmix_one - real(25,cvmix_r8)*zetastar)**(cvmix_one/real(3,cvmix_r8))
        end if
      end if
    else    ! not lStokesMOST

      if (lm) then
        if (zeta.ge.cvmix_zero) then
          ! Stable region
          compute_phi_inv = cvmix_one/(cvmix_one + real(5,cvmix_r8)*zeta)
        else if (zeta.ge.CVmix_kpp_params_in%zeta_m) then
          compute_phi_inv = (cvmix_one - real(16,cvmix_r8)*zeta)**0.25_cvmix_r8
        else
          compute_phi_inv = (CVmix_kpp_params_in%a_m -                          &
                            CVmix_kpp_params_in%c_m*zeta)**                     &
                            (cvmix_one/real(3,cvmix_r8))
        end if
      end if

      if (ls) then
        if (zeta.ge.cvmix_zero) then
          ! Stable region
          compute_phi_inv = cvmix_one/(cvmix_one + real(5,cvmix_r8)*zeta)
        else if (zeta.ge.CVmix_kpp_params_in%zeta_s) then
          compute_phi_inv = (cvmix_one - real(16,cvmix_r8)*zeta)**0.5_cvmix_r8
        else
          compute_phi_inv = (CVmix_kpp_params_in%a_s -                          &
                            CVmix_kpp_params_in%c_s*zeta)**                     &
                            (cvmix_one/real(3,cvmix_r8))
        end if
      end if
    end if

!EOC

  end function compute_phi_inv

!BOP

! !IROUTINE: compute_Stokes_chi
! !INTERFACE:

  function compute_Stokes_chi( xi, lchi_m, lchi_s)

! !DESCRIPTION:
!  Compute Stokes similarity function chi, of Stokes parameter xi= Ps/(PU+PS+PB)
!\\
!\\

! !INPUT PARAMETERS:
    real(cvmix_r8),              intent(in) :: xi
    logical, optional,           intent(in) :: lchi_m, lchi_s

! !OUTPUT PARAMETERS:
    real(cvmix_r8) :: compute_Stokes_chi

!EOP
!BOC

    logical :: lm, ls
    real(cvmix_r8)  :: chi, xi_tmp

    ! If not specifying lchi_m or lchi_s, routine will error out, but
    ! initializing result to 1 (No Stokes value) removes warning about
    ! possibly returning an un-initialized value
    compute_Stokes_chi = cvmix_one

  if (present(lchi_m)) then
      lm = lchi_m
  else
      lm = .false.
  end if

  if (present(lchi_s)) then
    ls = lchi_s
  else
    ls = .false.
  end if

  if (lm.eqv.ls) then
    print*, "ERROR: must compute chi_m or chi_s, can not compute both!"
    stop 1
  end if

  if (lm) then
    xi_tmp = MAX( cvmix_zero , MIN( xi , 0.73_cvmix_r8 ) )
    chi = cvmix_one - 1.671_cvmix_r8 * xi_tmp
    if ( xi_tmp .ge. 0.35_cvmix_r8) &
      chi = 1.03_cvmix_r8 + xi_tmp * ( 1.58_cvmix_r8*xi_tmp - 2.31_cvmix_r8 )
    compute_Stokes_chi = chi
  end if

  if (ls) then
    xi_tmp = MAX( cvmix_zero , MIN( xi , 0.89_cvmix_r8 ) )
    chi = cvmix_one - 1.594_cvmix_r8 * xi_tmp
    if ( xi_tmp .ge. 0.35_cvmix_r8) &
      chi = 0.78_cvmix_r8 + xi_tmp * ( 0.67_cvmix_r8*xi_tmp - 1.20_cvmix_r8 )
    compute_Stokes_chi = chi
  end if

!EOC

  end function compute_Stokes_chi

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

! !INPUT PARAMETERS:
    real(cvmix_r8), intent(in) :: GAT1  ! G(1)
    real(cvmix_r8), intent(in) :: DGAT1 ! G'(1)

! !OUTPUT PARAMETERS:
    real(cvmix_r8), dimension(4), intent(inout) :: coeffs

!EOP
!BOC

    coeffs(1) =  cvmix_zero
    coeffs(2) =  cvmix_one
    coeffs(3) =  real(3,cvmix_r8)*GAT1 - DGAT1 - real(2,cvmix_r8)
    coeffs(4) = -real(2,cvmix_r8)*GAT1 + DGAT1 + cvmix_one

!EOC

  end subroutine cvmix_kpp_compute_shape_function_coeffs

!BOP

! !IROUTINE: cvmix_compute_nu_at_OBL_depth_LMD94
! !INTERFACE:

  function cvmix_kpp_compute_nu_at_OBL_depth_LMD94(depths_cntr, layer_widths, &
                                                   diffs_iface, OBL_depth,    &
                                                   diff_2above, dnu_dz)

! !DESCRIPTION:
!  Interpolate to find $\nu$ at \verb|OBL_depth| from values at interfaces
!  above and below.
!\\
!\\

! !INPUT PARAMETERS:
    ! depths_cntr  = (/layer center containing OBL, layer center below/)
    ! diffs_iface  = diffusivity at interfaces of cell containing OBL
    ! layer_widths = (/width of layer containing OBL, width of layer below/)
    real(cvmix_r8), dimension(2), intent(in) :: depths_cntr, diffs_iface,     &
                                                layer_widths
    real(cvmix_r8),               intent(in) :: OBL_depth
    ! diffusivity at iface above the iface above OBL_depth (not needed if
    ! OBL is in top layer)
    real(cvmix_r8), optional,     intent(in) :: diff_2above

! !OUTPUT PARAMETERS:
    real(cvmix_r8), optional, intent(out) :: dnu_dz
    real(cvmix_r8) :: cvmix_kpp_compute_nu_at_OBL_depth_LMD94

!EOP
!BOC

    ! Local variables
    real(cvmix_r8), dimension(4) :: coeffs
    real(cvmix_r8) :: dnu_dz_above, dnu_dz_below, dnu_dz_local, wgt
    real(cvmix_r8) :: iface_depth

    ! (1) Compute derivatives of nu at layer centers (central difference)
    !     Sign convention: dnu/dz is positive if nu increases as you
    !                      move up in the column
    if (present(diff_2above)) then
      dnu_dz_above = (diff_2above-diffs_iface(1))/layer_widths(1)
    else
      ! Assume diffusivity goes to 0 at surface (z=0)
      dnu_dz_above = -diffs_iface(1)/layer_widths(1)
    end if
    dnu_dz_below = (diffs_iface(1)-diffs_iface(2))/layer_widths(2)
    ! Stability => require non-negative dnu_dz
    if (dnu_dz_above.lt.0.0_cvmix_r8) dnu_dz_above = 0.0_cvmix_r8
    if (dnu_dz_below.lt.0.0_cvmix_r8) dnu_dz_below = 0.0_cvmix_r8

    ! (2) Compute dnu/dz at OBL_depth by weighted average of values
    !     computed above (see LMD94, Eq. (D5) for details)
    iface_depth = depths_cntr(1) - 0.5_cvmix_r8*layer_widths(1)
    wgt = (-iface_depth-OBL_depth) / layer_widths(1)
    dnu_dz_local = wgt*dnu_dz_above + (cvmix_one-wgt)*dnu_dz_below

    ! (3) Linear interpolant: slope = value computed in (2) and the line goes
    !     through the point (iface_depth, diffs_iface(1))
    coeffs = cvmix_zero
    coeffs(1) = diffs_iface(1) - dnu_dz_local*iface_depth
    coeffs(2) = dnu_dz_local
    if (present(dnu_dz)) then
      dnu_dz = dnu_dz_local
    end if
    cvmix_kpp_compute_nu_at_OBL_depth_LMD94=cvmix_math_evaluate_cubic(coeffs, &
                                                                  -OBL_depth)
!EOC

  end function cvmix_kpp_compute_nu_at_OBL_depth_LMD94

!BOP

! !IROUTINE: cvmix_kpp_EFactor_model
! !INTERFACE:
  function cvmix_kpp_EFactor_model(u10, ustar, hbl, CVmix_params_in)

! !DESCRIPTION:
!  This function returns the enhancement factor, given the 10-meter
!  wind (m/s), friction velocity (m/s) and the boundary layer depth (m).
!\\
!\\

! !INPUT PARAMETERS:
    real(cvmix_r8),                 intent(in) :: u10    ! 10 meter wind (m/s)
    real(cvmix_r8),                 intent(in) :: ustar  ! water-side surface friction velocity (m/s)
    real(cvmix_r8),                 intent(in) :: hbl    ! boundary layer depth (m)
    type(cvmix_global_params_type), intent(in) :: CVmix_params_in

!EOP
!BOC

    real(cvmix_r8) :: us_sl, lasl_sqr_i
    real(cvmix_r8) :: cvmix_kpp_EFactor_model

    if (u10 .gt. cvmix_zero .and. ustar .gt. cvmix_zero) then
      ! surface layer averaged Stokes drift
      us_sl = cvmix_kpp_ustokes_SL_model(u10, hbl, CVmix_params_in)
      !
      ! LaSL^{-2}
      lasl_sqr_i = us_sl/ustar
      !
      ! enhancement factor (Li et al., 2016)
      cvmix_kpp_EFactor_model = sqrt(cvmix_one &
                 +cvmix_one/1.5_cvmix_r8**2*lasl_sqr_i &
                 +cvmix_one/5.4_cvmix_r8**4*lasl_sqr_i**2)
    else
      ! otherwise set to one
      cvmix_kpp_EFactor_model = cvmix_one
    end if

!EOC

  end function cvmix_kpp_EFactor_model

!BOP

! !IROUTINE: cvmix_kpp_ustokes_SL_model
! !INTERFACE:
  function cvmix_kpp_ustokes_SL_model(u10, hbl, CVmix_params_in)

! !DESCRIPTION:
!  This function returns the surface layer averaged Stokes drift, given
!  the 10-meter wind (m/s) and the boundary layer depth (m).
!\\
!\\

! !INPUT PARAMETERS:
    real(cvmix_r8),                 intent(in) :: u10    ! 10 meter wind (m/s)
    real(cvmix_r8),                 intent(in) :: hbl    ! boundary layer depth (m)
    type(cvmix_global_params_type), intent(in) :: CVmix_params_in

!EOP
!BOC

    real(cvmix_r8), parameter :: &
        ! ratio of U19.5 to U10 (Holthuijsen, 2007)
        u19p5_to_u10 = 1.075_cvmix_r8, &
        ! ratio of mean frequency to peak frequency for
        ! Pierson-Moskowitz spectrum (Webb, 2011)
        fm_to_fp = 1.296_cvmix_r8, &
        ! ratio of surface Stokes drift to U10
        us_to_u10 = 0.0162_cvmix_r8, &
        ! loss ratio of Stokes transport
        r_loss = 0.667_cvmix_r8

    real(cvmix_r8) :: us, hm0, fm, fp, vstokes, kphil, kstar
    real(cvmix_r8) :: z0, z0i, r1, r2, r3, r4, tmp
    real(cvmix_r8) :: cvmix_kpp_ustokes_SL_model

    if (u10 .gt. cvmix_zero) then
      ! surface Stokes drift
      us = us_to_u10*u10
      !
      ! significant wave height from Pierson-Moskowitz
      ! spectrum (Bouws, 1998)
      hm0 = 0.0246_cvmix_r8*u10**2
      !
      ! peak frequency (PM, Bouws, 1998)
      tmp = 2.0_cvmix_r8*cvmix_PI*u19p5_to_u10*u10
      fp = 0.877_cvmix_r8*CVmix_params_in%Gravity/tmp
      !
      ! mean frequency
      fm = fm_to_fp*fp
      !
      ! total Stokes transport (a factor r_loss is applied to account
      !  for the effect of directional spreading, multidirectional waves
      !  and the use of PM peak frequency and PM significant wave height
      !  on estimating the Stokes transport)
      vstokes = 0.125_cvmix_r8*cvmix_PI*r_loss*fm*hm0**2
      !
      ! the general peak wavenumber for Phillips' spectrum
      ! (Breivik et al., 2016) with correction of directional spreading
      kphil = 0.176_cvmix_r8*us/vstokes
      !
      ! surface layer averaged Stokes dirft with Stokes drift profile
      ! estimated from Phillips' spectrum (Breivik et al., 2016)
      ! the directional spreading effect from Webb and Fox-Kemper, 2015
      ! is also included
      kstar = kphil*2.56_cvmix_r8
      ! surface layer
      z0 = 0.2_cvmix_r8*abs(hbl)
      z0i = cvmix_one/z0
      ! term 1 to 4
      r1 = (0.151_cvmix_r8/kphil*z0i-0.84_cvmix_r8) &
            *(cvmix_one-exp(-2.0_cvmix_r8*kphil*z0))
      r2 = -(0.84_cvmix_r8+0.0591_cvmix_r8/kphil*z0i) &
             *sqrt(2.0_cvmix_r8*cvmix_PI*kphil*z0) &
             *erfc(sqrt(2.0_cvmix_r8*kphil*z0))
      r3 = (0.0632_cvmix_r8/kstar*z0i+0.125_cvmix_r8) &
            *(cvmix_one-exp(-2.0_cvmix_r8*kstar*z0))
      r4 = (0.125_cvmix_r8+0.0946_cvmix_r8/kstar*z0i) &
             *sqrt(2.0_cvmix_r8*cvmix_PI*kstar*z0) &
             *erfc(sqrt(2.0_cvmix_r8*kstar*z0))
      cvmix_kpp_ustokes_SL_model = us*(0.715_cvmix_r8+r1+r2+r3+r4)
    else
      cvmix_kpp_ustokes_SL_model = cvmix_zero
    end if

!EOC

  end function cvmix_kpp_ustokes_SL_model

!BOP

! !IROUTINE: cvmix_kpp_composite_shape
! !INTERFACE:

  function cvmix_kpp_composite_shape( sigma , Gat1)

!  !DESCRIPTION:
!   This function returns the value of the composite shape function for both
!   momentum and scalars at fractional depth sigma in the boundary layer.
!   This shape function is a cubic for sigma<sig\_m; and a quadratic below, as
!   fit to Fig. 6 of Large et al., 2020 (doi:10.1175/JPO-D-20-0308.1)
!   The subroutine also returns the derivative  dG / dsig
!\\
!\\

! !INPUT PARAMETERS:
    real(cvmix_r8),                    intent(in)  ::  sigma
    real(cvmix_r8), optional,          intent(in)  ::  Gat1

! !OUTPUT PARAMETERS:
    real(cvmix_r8) :: cvmix_kpp_composite_shape

!EOP
!BOC

    real(cvmix_r8)  :: a2Gsig, a3Gsig, sig_m, G_m, G_1, sig
    a2Gsig = -2.1637_cvmix_r8
    a3Gsig =  0.5831_cvmix_r8
    sig_m  =  0.35_cvmix_r8
    G_m    =  0.11_cvmix_r8     ! sig_m + sig_m * sig_m * (a2Gsig + a3Gsig * sig_m)
    if ( present(Gat1) ) then
      G_1  = MAX( cvmix_zero  ,  MIN( Gat1 , G_m ) )
    else
      G_1  =  cvmix_zero
    end if

    if (sigma .lt. sig_m)  then
      sig = MAX( sigma , cvmix_zero)
      cvmix_kpp_composite_shape = sig + sig * sig * (a2Gsig + a3Gsig * sig)
    else
      sig = MIN( sigma , cvmix_one )
      cvmix_kpp_composite_shape = G_1 + (G_m-G_1) * ((1.-sig) / (1.-sig_m))**2
    end if

!EOC

  end function cvmix_kpp_composite_shape

!BOP

! !IROUTINE: cvmix_kpp_composite_Gshape
! !INTERFACE:

  subroutine cvmix_kpp_composite_Gshape(sigma , Gat1, Gsig, dGdsig)

! !INPUT PARAMETERS:
    real(cvmix_r8),                    intent(in)  ::  sigma
    real(cvmix_r8),                    intent(in)  ::  Gat1

! !INPUT/OUTPUT PARAMETERS:
    real(cvmix_r8),                  intent(inout) ::  Gsig
    real(cvmix_r8),                  intent(inout) ::  dGdsig

!EOP
!BOC

    real(cvmix_r8)  :: a2Gsig, a3Gsig, bGsig, sig_m, G_m, G_1, sig

    a2Gsig = -2.1637_cvmix_r8
    a3Gsig =  0.5831_cvmix_r8
    sig_m  =  0.35_cvmix_r8
    G_m    =  0.11_cvmix_r8     ! sig_m + sig_m * sig_m * (a2Gsig + a3Gsig * sig_m)
    G_1    =  MAX( cvmix_zero  ,  MIN( Gat1 , G_m ) )

    if (sigma .lT. sig_m)  then
      sig    = MAX( sigma , cvmix_zero)
      Gsig   = sig + sig * sig * (a2Gsig + a3Gsig * sig)
      dGdsig = cvmix_one + sig * (2.0_cvmix_r8 * a2Gsig + 3.0_cvmix_r8 * a3Gsig * sig)
    else
      bGsig  =  (G_m-G_1) / (1.-sig_m)**2
      sig    = MIN( sigma , cvmix_one )
      Gsig   = G_1 + bGsig * (cvmix_one - sig) * (cvmix_one - sig)
      dGdsig = bGsig * 2.0_cvmix_r8 * (sig - cvmix_one)
    end if

!EOC

  end subroutine cvmix_kpp_composite_Gshape

!BOP

! !IROUTINE: cvmix_coeffs_bkgnd_wrap
! !INTERFACE:

  subroutine cvmix_kpp_compute_StokesXi (zi, zk, kSL, SLDepth,          &
           surf_buoy_force, surf_fric_vel, omega_w2x, uE, vE, uS, vS,   &
           uSbar, vSbar, uS_SLD, vS_SLD, uSbar_SLD, vSbar_SLD,          &
           StokesXI, CVmix_kpp_params_user)

! !DESCRIPTION:
!  Compute the Stokes similarity parameter, StokesXI, and Entrainment Rule, BEdE\_ER, from
!  surface layer integrated TKE production terms as parameterized in
!  Large et al., 2020 (doi:10.1175/JPO-D-20-0308.1)
!\\
!\\

! !INPUT PARAMETERS:
    real(cvmix_r8), dimension(:), intent(in) :: zi, zk          !< Cell interface and center heights < 0
    integer,        intent(in) :: kSL                           !< cell index of Surface Layer Depth
    real(cvmix_r8), intent(in) :: SLDepth                       !< Surface Layer Depth Integration limit
    real(cvmix_r8), intent(in) :: surf_buoy_force               !< Surface buoyancy flux forcing
    real(cvmix_r8), intent(in) :: surf_fric_vel, omega_w2x      !< Surface wind forcing from x-axis
    real(cvmix_r8), dimension(:), intent(in) :: uE, vE          !< Eulerian velocity at centers
    real(cvmix_r8), intent(in) :: uS_SLD, vS_SLD                !< Stokes drift at SLDepth
    real(cvmix_r8), intent(in) :: uSbar_SLD, vSbar_SLD          !< Average Stokes drift cell kSL to SLDepth

    type(cvmix_kpp_params_type), intent(in), optional, target :: CVmix_kpp_params_user

! !INPUT/OUTPUT PARAMETERS:
    real(cvmix_r8), dimension(:), intent(inout) :: uS, vS       !< Stokes drift at interfaces
    real(cvmix_r8), dimension(:), intent(inout) :: uSbar, vSbar !< Cell average Stokes drift
    real(cvmix_r8), intent(inout) :: StokesXI                   !< Stokes similarity parameter

!EOP
!BOC

    type(cvmix_kpp_params_type),  pointer :: CVmix_kpp_params_in
    real(cvmix_r8) :: PU, PS , PB                   ! surface layer TKE production terms
    real(cvmix_r8) :: uS_TMP, vS_TMP, uSbar_TMP, vSbar_TMP                       ! Temporary Store
    real(cvmix_r8) :: ustar, delH, delU, delV, omega_E2x, cosOmega, sinOmega
    real(cvmix_r8) :: BLDepth, TauMAG, TauCG, TauDG, taux0, tauy0, Stk0 , Pinc
    real(cvmix_r8) :: PBfact , CempCGm                                           ! Empirical constants
    real(cvmix_r8) :: dtop, tauEtop, tauxtop, tauytop                            ! Cell top values
    real(cvmix_r8) :: dbot, tauEbot, tauxbot, tauybot, sigbot, Gbot              ! Cell bottom values
    integer        :: ktmp                                                       ! vertical loop

    CVmix_kpp_params_in => CVmix_kpp_params_saved
    if (present(CVmix_kpp_params_user)) then
      CVmix_kpp_params_in => CVmix_kpp_params_user
    end if

    if ( CVmix_kpp_params_in%lStokesMOST ) then

    ! Move bottom of cell kSL up to Surface Layer Extent = SLDepth
    uS_TMP = uS(kSL+1)
    vS_TMP = vS(kSL+1)
    uSbar_TMP = uSbar(kSL)
    vSbar_TMP = vSbar(kSL)
    uS(kSL+1)   = uS_SLD
    vS(kSL+1)   = vS_SLD
    uSbar(kSL)= uSbar_SLD
    vSbar(kSL)= vSbar_SLD

    CempCGm=  3.5_cvmix_r8

    ustar   = MAX( surf_fric_vel , 1.e-4_cvmix_r8 )   ! > 0
    taux0   = ustar**2 * cos(omega_w2x)
    tauy0   = ustar**2 * sin(omega_w2x)
    Stk0    = sqrt( uS(1)**2 + vS(1)**2 )
    BLDepth = SLDepth / CVmix_kpp_params_in%surf_layer_ext

    ! Parameterized Buoyancy production of TKE
    PBfact =  0.110_cvmix_r8
    PB      = PBfact * MAX( -surf_buoy_force * BLdepth ,  cvmix_zero )

    ! Compute Both Shear Production Terms down from  Surface = initial top values
    PU      = 0.0
    PS      = 0.0
    dtop    = 0.0
    delU    = uE(1) - uE(2)
    delV    = vE(1) - vE(2)
    tauEtop = (taux0 * delU + tauy0 * delV ) / (zk(1) - zk(2) )
    tauxtop = taux0
    tauytop = tauy0

    do ktmp = 1, kSL
      ! SLdepth can be between cell interfaces kSL and kSL+1
      delH = min( max(cvmix_zero, SLdepth - dtop), (zi(ktmp) - zi(ktmp+1) ) )
      dbot = MIN( dtop + delH ,  SLdepth)
      sigbot = dbot / BLdepth
      Gbot     = cvmix_kpp_composite_shape(sigbot)
      TauMAG   = ustar * ustar * Gbot / sigbot
      delU     = uE(ktmp) - uE(ktmp+1)
      delV     = vE(ktmp) - vE(ktmp+1)
      Omega_E2x= atan2( delV  , delU )
      cosOmega = cos(Omega_E2x)
      sinOmega = sin(Omega_E2x)
      tauCG    = CempCGm * Gbot *  (taux0 * cosOmega - tauy0 * sinOmega)
      ! tauDG    = sqrt( TauMAG**2 - tauCG**2 ) ! G
      tauDG    = TauMAG                       ! E
      tauxbot  = tauDG * cosOmega  -  tauCG * sinOmega
      tauybot  = tauDG * sinOmega  +  tauCG * cosOmega
      tauEbot  = (tauxbot * delU + tauybot * delV) / (zk(ktmp) - zk(ktmp+1) )

      ! Increment Eulerian Shear Production
      Pinc     = 0.5_cvmix_r8 * (tauEbot + tauEtop) * delH
      PU       = PU + MAX( Pinc , cvmix_zero )

      ! Increment Stokes Shear Production
      Pinc   = tauxtop*uS(ktmp) - tauxbot*uS(ktmp+1) + tauytop*vS(ktmp) - tauybot*vS(ktmp+1)
      Pinc   = Pinc - (tauxtop-tauxbot) * uSbar(ktmp)  - (tauytop-tauybot) * vSbar(ktmp)
      PS     = PS +  MAX( Pinc , cvmix_zero )

      ! Bottom becomes next top
      dtop    = dbot
      tauxtop = tauxbot
      tauytop = tauybot
      tauEtop = tauEbot
    enddo

    ! Compute Stokes similarity parameter
    StokesXI  = PS / MAX( PU + PS + PB , 1.e-12_cvmix_r8 )


    ! Restore bottom of cell kSL at zi(kSL+1) with stored Stokes Drift ; ditto average over cell kSL
    uS(kSL+1)  = uS_TMP
    vS(kSL+1)  = vS_TMP
    uSbar(kSL) = uSbar_TMP
    vSbar(kSL) = vSbar_TMP

   else    ! not lStokesMOST
     StokesXI  = cvmix_zero
   end if

!EOC

  end subroutine cvmix_kpp_compute_StokesXi

end module cvmix_kpp
