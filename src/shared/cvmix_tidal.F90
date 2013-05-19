!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module cvmix_tidal

!BOP
!\newpage
! !MODULE: cvmix_tidal
!
! !DESCRIPTION:
!  This module contains routines to initialize the derived types needed for
!  tidal mixing (currently just the Simmons scheme) and to set the viscosity
!  and diffusivity coefficients accordingly.
!\\
!\\
!
! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !USES:

   use cvmix_kinds_and_types, only : cvmix_r8,                &
                                     cvmix_data_type,         &
                                     cvmix_tidal_params_type
!EOP

   implicit none
   private
   save

!BOP

! !PUBLIC MEMBER FUNCTIONS:

   public :: cvmix_init_tidal
   public :: cvmix_coeffs_tidal
!EOP

 contains

!BOP

! !IROUTINE: cvmix_init_tidal
! !INTERFACE:

  subroutine cvmix_init_tidal(CVmix_tidal_params, mix_scheme, units, &
                              efficiency, vertical_decay_scale,      &
                              max_coefficient, local_mixing_frac,    &
                              depth_cutoff)

! !DESCRIPTION:
!  Initialization routine for tidal mixing. There is currently just one
!  supported schemes - set \verb|mix_scheme = 'simmons'| to use the Simmons
!  mixing scheme.
!
! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*),         intent(in) :: mix_scheme
    character(len=*),         intent(in) :: units
    real(cvmix_r8), optional, intent(in) :: efficiency
    real(cvmix_r8), optional, intent(in) :: vertical_decay_scale
    real(cvmix_r8), optional, intent(in) :: max_coefficient
    real(cvmix_r8), optional, intent(in) :: local_mixing_frac
    real(cvmix_r8), optional, intent(in) :: depth_cutoff

! !OUTPUT PARAMETERS:
    type(cvmix_tidal_params_type), intent(inout) :: CVmix_tidal_params
!EOP
!BOC

    select case (trim(mix_scheme))
      case ('simmons','Simmons')
        CVmix_tidal_params%mix_scheme = trim(mix_scheme)

        ! Unitless parameters
        if (present(efficiency)) then
          CVmix_tidal_params%efficiency = efficiency
        else
          CVmix_tidal_params%efficiency = 0.2_cvmix_r8
        end if
        if (present(local_mixing_frac)) then
          CVmix_tidal_params%local_mixing_frac = local_mixing_frac
        else
          CVmix_tidal_params%local_mixing_frac = 1.0_cvmix_r8/3.0_cvmix_r8
        end if

        ! Parameters with units
        if (present(vertical_decay_scale)) then
          CVmix_tidal_params%vertical_decay_scale = vertical_decay_scale
        end if
        if (present(max_coefficient)) then
          CVmix_tidal_params%max_coefficient = max_coefficient
        end if
        if (present(depth_cutoff)) then
          CVmix_tidal_params%depth_cutoff = depth_cutoff
        else
          ! Default: no cutoff depth => 0 cm or 0 m
          CVmix_tidal_params%depth_cutoff = 0.0_cvmix_r8
        end if
        select case (trim(units))
          case ('mks')
            if (.not.present(vertical_decay_scale)) then
              CVmix_tidal_params%vertical_decay_scale = 500.0_cvmix_r8
            end if
            if (.not.present(max_coefficient)) then
              CVmix_tidal_params%max_coefficient = 50.0e-4_cvmix_r8
            end if

          case ('cgs')
            if (.not.present(vertical_decay_scale)) then
              CVmix_tidal_params%vertical_decay_scale = 500.0e2_cvmix_r8
            end if
            if (.not.present(max_coefficient)) then
              CVmix_tidal_params%max_coefficient = 50.0_cvmix_r8
            end if

          case DEFAULT
            print*, "ERROR: ", trim(units), " is not a valid choice for ",    &
                    "tidal mixing. Only 'mks' and 'cgs' are supported."
            stop 1

        end select

      case DEFAULT
        print*, "ERROR: ", trim(mix_scheme), " is not a valid choice for ", &
                "tidal mixing."
        stop 1

    end select

!EOC

  end subroutine cvmix_init_tidal

!***********************************************************************
!BOP
! !IROUTINE: cvmix_coeffs_tidal
! !INTERFACE:

  subroutine cvmix_coeffs_tidal(CVmix_vars, CVmix_tidal_params)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for tidal mixing
!  parameterizatiions.
!\\
!\\
!
! !USES:
!  only those used by entire module.

! !INPUT PARAMETERS:
    type(cvmix_tidal_params_type), intent(in) :: CVmix_tidal_params

! !INPUT/OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars

!EOP
!BOC

    ! Local variables
    real(cvmix_r8), allocatable, dimension(:) :: dep_fun ! Deposition Function
    integer :: nlev !, k
    real(cvmix_r8) :: surf_hgt, ocn_depth, decay_scale

    nlev      = CVmix_vars%nlev
    surf_hgt  = CVmix_vars%surf_hgt
    ocn_depth = CVmix_vars%ocn_depth
    decay_scale = CVmix_tidal_params%vertical_decay_scale
    allocate(dep_fun(nlev+1))
    ! Need to read / set z_iface before doing this loop!
!    do k=1,nlev+1
!      dep_fun(k) = exp(-CVmix_vars%z_iface(k)/decay_scale) !/ &
!         decay_scale*(exp(ocn_depth/decay_scale)-exp(-surf_hgt/decay_scale))
!    end do

    select case (trim(CVmix_tidal_params%mix_scheme))
      case ('simmons','Simmons')
          CVmix_vars%visc_iface = 0.0_cvmix_r8
          CVmix_vars%diff_iface = 0.0_cvmix_r8

      case DEFAULT
        ! Note: this error should be caught in cvmix_init_tidal
        print*, "ERROR: invalid choice for type of tidal mixing."
        stop 1

    end select

    deallocate(dep_fun)

!EOC
  end subroutine cvmix_coeffs_tidal

end module cvmix_tidal
