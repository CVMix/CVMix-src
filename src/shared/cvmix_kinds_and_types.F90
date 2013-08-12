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
  integer, parameter, public :: cvmix_r8     = selected_real_kind(15, 307), &
                                cvmix_strlen = 256

  ! Global parameters:
  ! The constant 1 is used repeatedly in PP and double-diff mixing.
  ! The value for pi is needed for Bryan-Lewis mixing.
  real(cvmix_r8), parameter, public :: one = 1.0_cvmix_r8
  real(cvmix_r8), parameter, public :: cvmix_PI = &
                                       3.14159265358979323846_cvmix_r8

! !PUBLIC TYPES:

  ! cvmix_data_type contains variables for time-dependent and column-specific
  ! mixing. Time-independent physical parameters should be stored in
  ! cvmix_global_params_type and *-mixing specific parameters should be
  ! stored in cvmix_*_params_type (found in the cvmix_* module).
  type, public :: cvmix_data_type
      integer        :: nlev = -1  ! Number of levels in column
                                   ! Setting default to -1 might be F95...

      ! Scalar quantities
      real(cvmix_r8) :: ocn_depth, & ! distance from sea level to ocean bottom
                                     ! (positive => below sea level)
                        OBL_depth, & ! distance from sea level to OBL bottom
                                     ! (positive => below sea level)
                        surf_hgt,  & ! sea surface height
                                     ! (positive => above sea level)
                        surf_fric, & ! turbulent friction velocity at surface
                        surf_buoy, & ! buoyancy forcing at surface
                        lat,       & ! latitude of column (degrees north)
                        lon,       & ! longitude of column (degrees east)
                        Coriolis,  & ! Coriolis parameter
                        kOBL_depth   ! index of cell containing OBL (fraction
                                     ! > .5 => below cell center)

      ! Values on interfaces
      ! For KPP, need to store non-local transport term
      ! nlev+1, 3
      real(cvmix_r8), dimension(:,:), pointer :: kpp_transport_iface => NULL()
      ! nlev+1, 2
      real(cvmix_r8), dimension(:,:), pointer :: diff_iface => NULL()
      ! nlev+1
      real(cvmix_r8), dimension(:),   pointer :: visc_iface => NULL()
      real(cvmix_r8), dimension(:),   pointer :: zw_iface   => NULL()
      real(cvmix_r8), dimension(:),   pointer :: dzw_iface  => NULL()
      real(cvmix_r8), dimension(:),   pointer :: Ri_iface   => NULL()
      ! For tidal mixing, we need the squared buoyancy frequency
      real(cvmix_r8), dimension(:),   pointer :: buoy_iface => NULL()

      ! Values at tracer points
      ! nlev
      real(cvmix_r8), dimension(:),   pointer :: dens     => NULL()
      real(cvmix_r8), dimension(:),   pointer :: dens_lwr => NULL()
      real(cvmix_r8), dimension(:),   pointer :: zt       => NULL()
      real(cvmix_r8), dimension(:),   pointer :: dzt      => NULL()
      real(cvmix_r8), dimension(:),   pointer :: Rib      => NULL()
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
      ! For densities, user must keep track of units (kg/m^3 vs g/cm^3)
      real(cvmix_r8)                 :: fw_rho    ! fresh water density
      real(cvmix_r8)                 :: sw_rho    ! salt water density
  end type cvmix_global_params_type

!EOP

end module cvmix_kinds_and_types

