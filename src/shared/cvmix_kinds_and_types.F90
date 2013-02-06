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
  ! The value for pi is needed for Bryan-Lewis mixing.
  real(cvmix_r8), parameter, public :: cvmix_PI = &
                      2.0_cvmix_r8*acos(0.0_cvmix_r8)

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
  end type cvmix_tidal_params_type

  ! cvmix_ddiff_params_type contains the necessary parameters for double
  ! diffusion mixing (currently just a place-holder variable)
  type, public :: cvmix_ddiff_params_type
      real(cvmix_r8) :: deleteme
  end type cvmix_ddiff_params_type

  ! cvmix_conv_params_type contains the necessary parameters for convective
  ! mixing.
  type, public :: cvmix_conv_params_type
      real(cvmix_r8)              :: convect_diff
      real(cvmix_r8)              :: convect_visc
  end type cvmix_conv_params_type
!EOP

end module cvmix_kinds_and_types

