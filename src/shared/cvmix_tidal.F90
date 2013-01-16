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

  subroutine cvmix_init_tidal(CVmix_tidal_params, mix_scheme)

! !DESCRIPTION:
!  Initialization routine for tidal mixing. There is currently just one
!  supported schemes - set \verb|mix_scheme = 'simmons'| to use the Simmons
!  mixing scheme.
!
! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: mix_scheme 

! !OUTPUT PARAMETERS:
    type(cvmix_tidal_params_type), intent(inout) :: CVmix_tidal_params
!EOP
!BOC

    select case (trim(mix_scheme))
      case ('simmons')
        CVmix_tidal_params%mix_scheme = "simmons"

      case DEFAULT
        print*, "ERROR: ", trim(mix_scheme), " is not a valid choice for ", &
                "tidal mixing."
        stop

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

    select case (trim(CVmix_tidal_params%mix_scheme))
      case ('simmons')
          CVmix_vars%visc_iface = 0.0_cvmix_r8
          CVmix_vars%diff_iface = 0.0_cvmix_r8

      case DEFAULT
        ! Note: this error should be caught in cvmix_init_tidal
        print*, "ERROR: invalid choice for type of tidal mixing."
        stop

    end select

!EOC
  end subroutine cvmix_coeffs_tidal

end module cvmix_tidal
