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

  ! Parameters to allow CVMix to store integers instead of strings
  integer, parameter, public :: CVMIX_OVERWRITE_OLD_VAL    = 1
  integer, parameter, public :: CVMIX_SUM_OLD_AND_NEW_VALS = 2
  integer, parameter, public :: CVMIX_MAX_OLD_AND_NEW_VALS = 3

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
                                     ! units: m
                        OBL_depth, & ! distance from sea level to OBL bottom
                                     ! (positive => below sea level)
                                     ! units: m
                        surf_hgt,  & ! sea surface height
                                     ! (positive => above sea level)
                                     ! units: m
                        surf_fric, & ! turbulent friction velocity at surface
                                     ! units: m/s
                        surf_buoy, & ! buoyancy forcing at surface
                                     ! units: m^2 s^-3
                        lat,       & ! latitude of column (degrees north)
                                     ! units: can be degrees or radians (there
                                     !        are no internal computation based
                                     !        on this term)
                        lon,       & ! longitude of column (degrees east)
                                     ! units: can be degrees or radians (there
                                     !        are no internal computation based
                                     !        on this term)
                        Coriolis,  & ! Coriolis parameter
                                     ! units: s^-1
                        kOBL_depth   ! index of cell containing OBL (fraction
                                     ! > .5 => below cell center)
                                     ! units: unitless

      ! Values on interfaces
      ! For KPP, need to store non-local transport term
      ! (:,1) = temperature tracer
      ! (:,2) = salinity / all non-temperature tracers
      ! (:,3) and (:,4) = momentum terms (x- and y-, respectively). Note that
      !                   currently both momentum terms are 0 everywhere
      ! Note that kpp_transport_iface is the value of K_x*gamma_x/flux_x: in
      ! other words, the user must multiply this value by either the freshwater
      ! flux or the penetrative shortwave heat flux to come the values in Eqs.
      ! (7.128) and (7.129) of the CVMix manual.
      ! nlev+1, 4
      real(cvmix_r8), dimension(:,:), pointer :: kpp_transport_iface => NULL()
                                             ! units: unitless (see note above)
      ! nlev+1, 2
      ! diffusivity coefficients at interfaces (2 columns needed for double diff)
      real(cvmix_r8), dimension(:,:), pointer :: diff_iface => NULL()
                                              ! units: m^2/s
      ! nlev+1
      ! viscosity (momentum diffusivity) coefficients at interfaces
      real(cvmix_r8), dimension(:),   pointer :: visc_iface => NULL()
                                              ! units: m^2/s
      ! height of interfaces in column (positive up => most are negative)
      real(cvmix_r8), dimension(:),   pointer :: zw_iface   => NULL()
                                              ! units: m
      ! distance between neighboring cell centers (first value is top of ocean
      ! to middle of first cell, last value is middle of last cell to ocean
      ! bottom)
      real(cvmix_r8), dimension(:),   pointer :: dzw_iface  => NULL()
                                              ! units: m
      ! shear Richardson number at column interfaces
      real(cvmix_r8), dimension(:),   pointer :: Ri_iface   => NULL()
                                              ! units: unitless
      ! For tidal mixing, we need the squared buoyancy frequency
      real(cvmix_r8), dimension(:),   pointer :: buoy_iface => NULL()
                                              ! units: s^-2

      ! Values at tracer points
      ! nlev
      ! Two density values are stored: the actual density of water (dens) and
      ! the density of water after adiabatic displacement to the level below
      ! where the water actually is (dens_lwr)
      real(cvmix_r8), dimension(:),   pointer :: dens     => NULL()
                                              ! units: kg m^-3 
      real(cvmix_r8), dimension(:),   pointer :: dens_lwr => NULL()
                                              ! units: kg m^-3 
      ! height of cell centers in column (positive up => most are negative)
      real(cvmix_r8), dimension(:),   pointer :: zt       => NULL()
                                              ! units: m
      ! level thicknesses (positive semi-definite)
      real(cvmix_r8), dimension(:),   pointer :: dzt      => NULL()
                                              ! units: m
      ! bulk Richardson number
      real(cvmix_r8), dimension(:),   pointer :: Rib      => NULL()
                                              ! units: unitless
      ! For double diffusion mixing, we need to calculate the stratification
      ! parameter R_rho. Since the denominator of this ratio may be zero,
      ! we store the numerator and denominator separately and make sure the
      ! denominator is non-zero before performing the division.
      real(cvmix_r8), dimension(:),   pointer :: strat_param_num   => NULL()
                                              ! units: unitless
      real(cvmix_r8), dimension(:),   pointer :: strat_param_denom => NULL()
                                              ! units: unitless
      ! For KPP we need buoyancy (as opposed to buoyancy frequency) and
      ! velocity (in both x direction and y direction)
      real(cvmix_r8), dimension(:),   pointer :: buoyancy          => NULL()
                                              ! units: m/(s^2)
      real(cvmix_r8), dimension(:),   pointer :: Vx                => NULL()
                                              ! units: m/s
      real(cvmix_r8), dimension(:),   pointer :: Vy                => NULL()
                                              ! units: m/s
  end type cvmix_data_type

  ! cvmix_global_params_type contains global parameters used by multiple
  ! mixing methods.
  type, public :: cvmix_global_params_type
      integer                        :: max_nlev  ! maximum number of levels
                                     ! units: unitless
      real(cvmix_r8)                 :: prandtl   ! Prandtl number
                                     ! units: unitless
      ! For densities, user must keep track of units (kg/m^3 vs g/cm^3)
      real(cvmix_r8)                 :: fw_rho    ! fresh water density
                                     ! units: kg m^-3 
      real(cvmix_r8)                 :: sw_rho    ! salt water density
                                     ! units: kg m^-3 
  end type cvmix_global_params_type

!EOP

end module cvmix_kinds_and_types

