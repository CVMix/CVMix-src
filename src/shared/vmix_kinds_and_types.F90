module vmix_kinds_and_types

!BOP
!
! !MODULE:  vmix_kinds_and_types
!
! !DESCRIPTION:
!  This module contains the declarations for all required vertical mixing
!  data types. It also contains several global parameters used by the vmix
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
  ! The vmix package uses double precision for floating point computations.
  integer, parameter, public :: vmix_r8     = selected_real_kind(13), &
                                vmix_strlen = 256

  ! Global parameters:
  ! The value for pi is needed for Bryan-Lewis mixing.
  real(kind=vmix_r8), parameter, public :: vmix_PI = &
                      2.0_vmix_r8*acos(0.0_vmix_r8)

! !PUBLIC TYPES:

  ! vmix_input_type contains every possible necessary input field for all
  ! supported types of mixing.
  type, public :: vmix_data_type
      integer :: nlev = -1 ! Number of levels in column
                           ! Setting default to -1 might be F95...

      ! Values on interfaces
      real(vmix_r8), dimension(:),   pointer :: visc_iface  ! nlev+1
      real(vmix_r8), dimension(:,:), pointer :: diff_iface  ! nlev+1, 2
      real(vmix_r8), dimension(:),   pointer :: z_iface     ! nlev+1
      real(vmix_r8), dimension(:),   pointer :: dw_iface    ! nlev+1

      ! Values at tracer points
      real(vmix_r8), dimension(:),   pointer :: dens        ! nlev
      real(vmix_r8), dimension(:),   pointer :: dens_lwr    ! nlev
      real(vmix_r8), dimension(:),   pointer :: z           ! nlev
      real(vmix_r8), dimension(:),   pointer :: dz          ! nlev
  end type vmix_data_type

  ! vmix_global_params_type contains global parameters used by multiple
  ! mixing methods.
  type, public :: vmix_global_params_type
      integer                       :: max_nlev  ! maximum number of levels
      real(vmix_r8)                 :: prandtl   ! Prandtl number
  end type vmix_global_params_type

  ! vmix_bkgnd_params_type contains the necessary parameters for background
  ! mixing. Background mixing fields can vary from level to level as well as
  ! over latitude and longitude.
  type, public :: vmix_bkgnd_params_type
      real(vmix_r8), allocatable :: static_visc(:,:) ! ncol, nlev+1
      real(vmix_r8), allocatable :: static_diff(:,:) ! ncol, nlev+1

      ! Note: need to include some logic to avoid excessive memory use
      !       when static_visc and static_diff are constant or 1-D
      logical :: lvary_vertical   ! True => second dim not 1
      logical :: lvary_horizontal ! True => first dim not 1
  end type vmix_bkgnd_params_type

  ! vmix_conv_params_type contains the necessary parameters for convective
  ! mixing.
  type, public :: vmix_conv_params_type
      real(vmix_r8)              :: convect_diff
      real(vmix_r8)              :: convect_visc
  end type vmix_conv_params_type

!EOP

end module vmix_kinds_and_types

