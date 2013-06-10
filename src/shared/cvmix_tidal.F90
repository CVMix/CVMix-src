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
                                     cvmix_strlen,             &
                                     cvmix_global_params_type
!EOP

   implicit none
   private
   save

!BOP

! !PUBLIC MEMBER FUNCTIONS:

   public :: cvmix_init_tidal
   public :: cvmix_compute_vert_dep
   public :: cvmix_coeffs_tidal
   public :: cvmix_put_tidal
   public :: cvmix_get_tidal_real
   public :: cvmix_get_tidal_str

   interface cvmix_put_tidal
     module procedure cvmix_put_tidal_real
     module procedure cvmix_put_tidal_str
   end interface cvmix_put_tidal

! !PUBLIC TYPES:

   ! cvmix_tidal_params_type contains the necessary parameters for tidal mixing
   ! (currently just Simmons)
   type, public :: cvmix_tidal_params_type
      private
      character(len=cvmix_strlen) :: mix_scheme
      real(cvmix_r8)              :: efficiency
      real(cvmix_r8)              :: vertical_decay_scale
      real(cvmix_r8)              :: max_coefficient
      real(cvmix_r8)              :: local_mixing_frac
      real(cvmix_r8)              :: depth_cutoff
   end type cvmix_tidal_params_type
!EOP

 contains

!BOP

! !IROUTINE: cvmix_init_tidal
! !INTERFACE:

  subroutine cvmix_init_tidal(CVmix_tidal_params, mix_scheme, units,          &
                              efficiency, vertical_decay_scale,               &
                              max_coefficient, local_mixing_frac, depth_cutoff)

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

!BOP

! !IROUTINE: cvmix_coeffs_tidal
! !INTERFACE:

  subroutine cvmix_coeffs_tidal(CVmix_vars, CVmix_tidal_params, CVmix_params, &
                                energy_flux)

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
    real(cvmix_r8),                 intent(in) :: energy_flux

! !INPUT/OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars

!EOP
!BOC

    ! Local variables
    integer        :: nlev, k
    real(cvmix_r8) :: coef, rho, buoy, z_cut
    real(cvmix_r8), allocatable, dimension(:) :: vert_dep

    nlev = CVmix_vars%nlev
    rho  = CVmix_params%fw_rho

    select case (trim(CVmix_tidal_params%mix_scheme))
      case ('simmons','Simmons')
        allocate(vert_dep(nlev+1))
        vert_dep = cvmix_compute_vert_dep(CVmix_vars, CVmix_tidal_params)
        coef = CVmix_tidal_params%local_mixing_frac * &
               CVmix_tidal_params%efficiency *        &
               energy_flux
        CVmix_vars%diff_iface = 0.0_cvmix_r8
        if (CVmix_vars%ocn_depth.ge.CVmix_tidal_params%depth_cutoff) then
          do k=1, nlev+1
            buoy = CVmix_vars%buoy_iface(k)
            z_cut = CVmix_tidal_params%depth_cutoff
            if (buoy.gt.0.0_cvmix_r8) &
              CVmix_vars%diff_iface(k,1) = coef*vert_dep(k)/(rho*buoy)
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

!BOP

! !IROUTINE: cvmix_compute_vert_dep
! !INTERFACE:

  function cvmix_compute_vert_dep(CVmix_vars, CVmix_tidal_params)

! !DESCRIPTION:
!  Computes the vertical deposition function needed for Simmons et al tidal
!  mixing.
!\\
!\\
!
! !USES:
!  only those used by entire module.

! !INPUT PARAMETERS:
    type(cvmix_tidal_params_type), intent(in) :: CVmix_tidal_params
    type(cvmix_data_type),         intent(in) :: CVmix_vars

! !OUTPUT PARAMETERS:
    real(cvmix_r8), dimension(CVMix_vars%nlev+1) :: cvmix_compute_vert_dep

!EOP
!BOC

    ! Local variables
    real(cvmix_r8) :: tot_area, num, thick
    integer        :: k, nlev

    nlev = CVmix_vars%nlev

    ! Compute vertical deposition
    tot_area = 0.0_cvmix_r8
    cvmix_compute_vert_dep(1) = 0.0_cvmix_r8
    cvmix_compute_vert_dep(nlev+1) = 0.0_cvmix_r8
    do k=2,nlev
      num = -CVmix_vars%zw_iface(k)/CVmix_tidal_params%vertical_decay_scale
      ! Simmons vertical deposition
      ! Note that it is getting normalized (divide through by tot_area)
      ! So multiplicative constants that are independent of z are omitted
      cvmix_compute_vert_dep(k) = exp(num)

      ! Compute integral of vert_dep via trapezoid rule
      ! (looks like midpoint rule, but vert_dep = 0 at z=0 and z=-ocn_depth)
      thick = CVmix_vars%zt(k-1) - CVmix_vars%zt(k)
      tot_area = tot_area + cvmix_compute_vert_dep(k)*thick
    end do
    ! Normalize vert_dep (need integral = 1.0D0)
    cvmix_compute_vert_dep = cvmix_compute_vert_dep/tot_area

!EOC

  end function cvmix_compute_vert_dep

!BOP

! !IROUTINE: cvmix_put_tidal_real
! !INTERFACE:

  subroutine cvmix_put_tidal_real(CVmix_tidal_params, varname, val)

! !DESCRIPTION:
!  Write a real value into a cvmix\_tidal\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    real(cvmix_r8),   intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_tidal_params_type), intent(inout) :: CVmix_tidal_params
!EOP
!BOC

    select case (trim(varname))
      case ('efficiency')
        CVmix_tidal_params%efficiency = val
      case ('vertical_decay_scale')
        CVmix_tidal_params%vertical_decay_scale = val
      case ('max_coefficient')
        CVmix_tidal_params%max_coefficient = val
      case ('local_mixing_frac')
        CVmix_tidal_params%local_mixing_frac = val
      case ('depth_cutoff')
        CVmix_tidal_params%depth_cutoff = val
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end subroutine cvmix_put_tidal_real

!BOP

! !IROUTINE: cvmix_put_tidal_str
! !INTERFACE:

  subroutine cvmix_put_tidal_str(CVmix_tidal_params, varname, val)

! !DESCRIPTION:
!  Write a string into a cvmix\_tidal\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_tidal_params_type), intent(inout) :: CVmix_tidal_params
!EOP
!BOC

    select case (trim(varname))
      case ('mix_scheme')
        CVmix_tidal_params%mix_scheme = val
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end subroutine cvmix_put_tidal_str

!BOP

! !IROUTINE: cvmix_get_tidal_real
! !INTERFACE:

  function cvmix_get_tidal_real(CVmix_tidal_params, varname)

! !DESCRIPTION:
!  Returns the real value of a cvmix\_tidal\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    type(cvmix_tidal_params_type), intent(in) :: CVmix_tidal_params

! !OUTPUT PARAMETERS:
    real(cvmix_r8) :: cvmix_get_tidal_real
!EOP
!BOC

    cvmix_get_tidal_real = 0.0_cvmix_r8
    select case (trim(varname))
      case ('efficiency')
        cvmix_get_tidal_real = CVmix_tidal_params%efficiency
      case ('vertical_decay_scale')
        cvmix_get_tidal_real = CVmix_tidal_params%vertical_decay_scale
      case ('max_coefficient')
        cvmix_get_tidal_real = CVmix_tidal_params%max_coefficient
      case ('local_mixing_frac')
        cvmix_get_tidal_real = CVmix_tidal_params%local_mixing_frac
      case ('depth_cutoff')
        cvmix_get_tidal_real = CVmix_tidal_params%depth_cutoff
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end function cvmix_get_tidal_real

!BOP

! !IROUTINE: cvmix_get_tidal_str
! !INTERFACE:

  function cvmix_get_tidal_str(CVmix_tidal_params, varname)

! !DESCRIPTION:
!  Returns the string value of a cvmix\_tidal\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    type(cvmix_tidal_params_type), intent(inout) :: CVmix_tidal_params

! !OUTPUT PARAMETERS:
    character(len=cvmix_strlen) :: cvmix_get_tidal_str
!EOP
!BOC

    select case (trim(varname))
      case ('mix_scheme')
        cvmix_get_tidal_str = CVmix_tidal_params%mix_scheme
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end function cvmix_get_tidal_str

end module cvmix_tidal
