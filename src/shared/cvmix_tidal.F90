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

   use cvmix_kinds_and_types, only : cvmix_r8,                 &
                                     cvmix_data_type,          &
                                     cvmix_tidal_params_type,  &
                                     cvmix_global_params_type
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

  subroutine cvmix_init_tidal(CVmix_tidal_params, CVmix_vars, mix_scheme, &
                              units, energy_flux, efficiency,             &
                              vertical_decay_scale, max_coefficient,      &
                              local_mixing_frac, depth_cutoff)

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
    real(cvmix_r8),           intent(in) :: energy_flux
    real(cvmix_r8), optional, intent(in) :: efficiency
    real(cvmix_r8), optional, intent(in) :: vertical_decay_scale
    real(cvmix_r8), optional, intent(in) :: max_coefficient
    real(cvmix_r8), optional, intent(in) :: local_mixing_frac
    real(cvmix_r8), optional, intent(in) :: depth_cutoff

! !OUTPUT PARAMETERS:
    type(cvmix_tidal_params_type), intent(inout) :: CVmix_tidal_params
    type(cvmix_data_type),         intent(inout) :: CVmix_vars
!EOP
!BOC

    ! Local variables
    real(cvmix_r8) :: tot_area, num, thick!, denom1, denom2, tmp
    integer        :: k, nlev

    nlev = CVmix_vars%nlev
    ! energy flux must be passed in (no default value exists)
    ! user is responsible for keeping track of units
    CVmix_tidal_params%energy_flux = energy_flux

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

    ! Compute vertical deposition function
    if(.not.associated(CVmix_vars%vert_dep)) then
      allocate(CVmix_vars%vert_dep(nlev+1))
    else
      if (size(CVmix_vars%vert_dep).ne.nlev+1) then
        write(*,"(A,1X,A,1X,I0)") "ERROR: vertical deposition function must", &
                                  "be array of length", nlev+1
        stop 1
      end if
    end if

    tot_area = 0.0_cvmix_r8
    CVmix_vars%vert_dep(1) = 0.0_cvmix_r8
    CVmix_vars%vert_dep(nlev+1) = 0.0_cvmix_r8
    do k=2,nlev
      num = -CVmix_vars%z_iface(k)/CVmix_tidal_params%vertical_decay_scale
      ! Simmons vertical deposition
      ! Note that it is getting normalized (divide through by tot_area)
      ! So multiplicative constants that are independent of z are omitted
      CVmix_vars%vert_dep(k) = exp(num)

      ! Compute integral of vert_dep via trapezoid rule
      ! (looks like midpoint rule, but vert_dep = 0 at z=0 and z=-ocn_depth)
      thick = CVmix_vars%z(k-1) - CVmix_vars%z(k)
      tot_area = tot_area + CVmix_vars%vert_dep(k)*thick
    end do
    ! Normalize vert_dep (need integral = 1.0D0)
    CVmix_vars%vert_dep = CVmix_vars%vert_dep/tot_area

!EOC

  end subroutine cvmix_init_tidal

!***********************************************************************
!BOP
! !IROUTINE: cvmix_coeffs_tidal
! !INTERFACE:

  subroutine cvmix_coeffs_tidal(CVmix_vars, CVmix_tidal_params, CVmix_params)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for tidal mixing
!  parameterizatiions.
!\\
!\\
!
! !USES:
!  only those used by entire module.

! !INPUT PARAMETERS:
    type(cvmix_tidal_params_type),  intent(in) :: CVmix_tidal_params
    type(cvmix_global_params_type), intent(in) :: CVmix_params

! !INPUT/OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars

!EOP
!BOC

    ! Local variables
    integer        :: nlev, k
    real(cvmix_r8) :: coef, rho, buoy, z_cut

    nlev = CVmix_vars%nlev
    rho  = CVmix_params%fw_rho

    select case (trim(CVmix_tidal_params%mix_scheme))
      case ('simmons','Simmons')
          coef = CVmix_tidal_params%local_mixing_frac*CVmix_tidal_params%efficiency*CVmix_tidal_params%energy_flux
          CVmix_vars%diff_iface = 0.0_cvmix_r8
          if (CVmix_vars%ocn_depth.ge.CVmix_tidal_params%depth_cutoff) then
            do k=1, nlev+1
              buoy = CVmix_vars%buoy(k)
              z_cut = CVmix_tidal_params%depth_cutoff
              if (buoy.gt.0.0_cvmix_r8) &
                CVmix_vars%diff_iface(k,1) = coef*CVmix_vars%vert_dep(k)/(rho*buoy)
              if (CVmix_vars%diff_iface(k,1).gt.CVmix_tidal_params%max_coefficient) &
                CVmix_vars%diff_iface(k,1) = CVmix_tidal_params%max_coefficient
            end do
          end if

      case DEFAULT
        ! Note: this error should be caught in cvmix_init_tidal
        print*, "ERROR: invalid choice for type of tidal mixing."
        stop 1

    end select

!EOC
  end subroutine cvmix_coeffs_tidal

end module cvmix_tidal
