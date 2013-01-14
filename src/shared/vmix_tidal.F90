!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module vmix_tidal

!BOP
!\newpage
! !MODULE: vmix_tidal
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

   use vmix_kinds_and_types, only : vmix_r8,                  &
                                    vmix_data_type,           &
                                    vmix_bkgnd_params_type,   &
                                    vmix_tidal_params_type
!EOP

   implicit none
   private
   save

!BOP

! !PUBLIC MEMBER FUNCTIONS:

   public :: vmix_init_tidal
   public :: vmix_coeffs_tidal
!EOP

 contains

!BOP

! !IROUTINE: vmix_init_tidal
! !INTERFACE:

  subroutine vmix_init_tidal(Vmix_tidal_params, mix_scheme)

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
    type(vmix_tidal_params_type), intent(inout) :: Vmix_tidal_params
!EOP
!BOC

    select case (trim(mix_scheme))
      case ('simmons')
        Vmix_tidal_params%mix_scheme = "simmons"

      case DEFAULT
        print*, "ERROR: ", trim(mix_scheme), " is not a valid choice for ", &
                "tidal mixing."
        stop

    end select

!EOC

  end subroutine vmix_init_tidal

!***********************************************************************
!BOP
! !IROUTINE: vmix_coeffs_tidal
! !INTERFACE:

  subroutine vmix_coeffs_tidal(Vmix_vars, Vmix_tidal_params)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for tidal mixing
!  parameterizatiions.
!\\
!\\
!
! !USES:
!  only those used by entire module.

! !INPUT PARAMETERS:
    type(vmix_tidal_params_type), intent(in) :: Vmix_tidal_params

! !INPUT/OUTPUT PARAMETERS:
    type(vmix_data_type), intent(inout) :: Vmix_vars
!EOP
!BOC

    select case (trim(Vmix_tidal_params%mix_scheme))
      case ('simmons')
          Vmix_vars%visc_iface = 0.0_vmix_r8
          Vmix_vars%diff_iface = 0.0_vmix_r8

      case DEFAULT
        ! Note: this error should be caught in vmix_init_tidal
        print*, "ERROR: invalid choice for type of tidal mixing."
        stop

    end select

!EOC
  end subroutine vmix_coeffs_tidal

end module vmix_tidal
