module cvmix_kinds_and_types

!BOP
!\newpage
! !MODULE:  cvmix_kinds_and_types
!
! !DESCRIPTION:
!  This module contains the declarations for all required vertical mixing
!  data types. It also contains several global parameters used by the cvmix
!  package, such as kind numbers and string lengths.
!  \\
!  \\

! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !USES:
!  uses no other modules
!EOP

  implicit none
  private
  save

!BOP

! !DEFINED PARAMETERS:

  ! Kind Types:
  ! The cvmix package uses double precision for floating point computations.
  integer, parameter, public :: cvmix_r8     = selected_real_kind(13), &
                                cvmix_strlen = 256

  ! Global parameters:
  ! The constant 1 is used repeatedly in PP and double-diff mixing.
  ! The value for pi is needed for Bryan-Lewis mixing.
  real(cvmix_r8), parameter, public :: one = 1.0_cvmix_r8
  real(cvmix_r8), parameter, public :: cvmix_PI = &
                                       3.14159265358979323846_cvmix_r8

! !PUBLIC TYPES:

  ! cvmix_input_type contains every possible necessary input field for all
  ! supported types of mixing.
  type, public :: cvmix_data_type
      integer :: nlev = -1 ! Number of levels in column
                           ! Setting default to -1 might be F95...

      ! Values on interfaces
      ! nlev+1, 2
      real(cvmix_r8), dimension(:,:), pointer :: diff_iface => NULL()
      ! nlev+1
      real(cvmix_r8), dimension(:),   pointer :: visc_iface => NULL()
      real(cvmix_r8), dimension(:),   pointer :: z_iface    => NULL()
      real(cvmix_r8), dimension(:),   pointer :: dw_iface   => NULL()
      real(cvmix_r8), dimension(:),   pointer :: Ri_iface => NULL()

      ! Values at tracer points
      ! nlev
      real(cvmix_r8), dimension(:),   pointer :: dens     => NULL()
      real(cvmix_r8), dimension(:),   pointer :: dens_lwr => NULL()
      real(cvmix_r8), dimension(:),   pointer :: z        => NULL()
      real(cvmix_r8), dimension(:),   pointer :: dz       => NULL()
      ! For double diffusion mixing, we need to calculate the stratification
      ! parameter R_rho. Since the denominator of this ratio may be zero,
      ! we store the numerator and denominator separately and make sure the
      ! denominator is non-zero before performing the division.
      real(cvmix_r8), dimension(:),   pointer :: strat_param_num   => NULL()
      real(cvmix_r8), dimension(:),   pointer :: strat_param_denom => NULL()
  end type cvmix_data_type

  ! cvmix_global_params_type contains global parameters used by multiple
  ! mixing methods.
  type, public :: cvmix_global_params_type
      integer                        :: max_nlev  ! maximum number of levels
      real(cvmix_r8)                 :: prandtl   ! Prandtl number
  end type cvmix_global_params_type

  ! cvmix_bkgnd_params_type contains the necessary parameters for background
  ! mixing. Background mixing fields can vary from level to level as well as
  ! over latitude and longitude.
  type, public :: cvmix_bkgnd_params_type
      real(cvmix_r8), allocatable :: static_visc(:,:) ! ncol, nlev+1
      real(cvmix_r8), allocatable :: static_diff(:,:) ! ncol, nlev+1

      ! Note: need to include some logic to avoid excessive memory use
      !       when static_visc and static_diff are constant or 1-D
      logical                     :: lvary_vertical   ! True => second dim not 1
      logical                     :: lvary_horizontal ! True => first dim not 1
  end type cvmix_bkgnd_params_type

  ! cvmix_shear_params_type contains the necessary parameters for shear mixing
  ! (currently Pacanowski-Philander or Large et al)
  type, public :: cvmix_shear_params_type
      character(len=cvmix_strlen) :: mix_scheme
      real(cvmix_r8)              :: PP_nu_zero
      real(cvmix_r8)              :: PP_alpha
      real(cvmix_r8)              :: PP_exp
      real(cvmix_r8)              :: KPP_nu_zero
      real(cvmix_r8)              :: KPP_Ri_zero
      real(cvmix_r8)              :: KPP_exp
  end type cvmix_shear_params_type

  ! cvmix_tidal_params_type contains the necessary parameters for shear mixing
  ! (currently just Simmons)
  type, public :: cvmix_tidal_params_type
      character(len=cvmix_strlen) :: mix_scheme
      real(cvmix_r8)              :: efficiency
      real(cvmix_r8)              :: vertical_decay_scale
      real(cvmix_r8)              :: max_coefficient
      real(cvmix_r8)              :: local_mixing_frac
      real(cvmix_r8)              :: depth_cutoff
  end type cvmix_tidal_params_type

  ! cvmix_ddiff_params_type contains the necessary parameters for double
  ! diffusion mixing
  type, public :: cvmix_ddiff_params_type
      real(cvmix_r8)              :: strat_param_max
      real(cvmix_r8)              :: kappa_ddiff_t
      real(cvmix_r8)              :: kappa_ddiff_s
      real(cvmix_r8)              :: ddiff_exp1
      real(cvmix_r8)              :: ddiff_exp2
      real(cvmix_r8)              :: kappa_ddiff_param1
      real(cvmix_r8)              :: kappa_ddiff_param2
      real(cvmix_r8)              :: kappa_ddiff_param3
      real(cvmix_r8)              :: mol_diff
  end type cvmix_ddiff_params_type

  ! cvmix_conv_params_type contains the necessary parameters for convective
  ! mixing.
  type, public :: cvmix_conv_params_type
      real(cvmix_r8)              :: convect_diff
      real(cvmix_r8)              :: convect_visc
  end type cvmix_conv_params_type
!EOP

end module cvmix_kinds_and_types

