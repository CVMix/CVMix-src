!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module vmix_ddiff

!BOP
!\newpage
! !MODULE: vmix_ddiff
!
! !DESCRIPTION:
!  This module contains routines to initialize the derived types needed for
!  double diffusion mixing and to set the viscosity and diffusivity
!  coefficients accordingly.
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
                                    vmix_ddiff_params_type
!EOP

   implicit none
   private
   save

!BOP

! !PUBLIC MEMBER FUNCTIONS:

   public :: vmix_init_ddiff
   public :: vmix_coeffs_ddiff
!EOP

 contains

!BOP

! !IROUTINE: vmix_init_ddiff
! !INTERFACE:

  subroutine vmix_init_ddiff(Vmix_ddiff_params)

! !DESCRIPTION:
!  Initialization routine for double diffusion mixing.
!
! !USES:
!  Only those used by entire module.

! !OUTPUT PARAMETERS:
    type(vmix_ddiff_params_type), intent(inout) :: Vmix_ddiff_params
!EOP
!BOC

    Vmix_ddiff_params%deleteme = 0.0_vmix_r8
!EOC

  end subroutine vmix_init_ddiff

!***********************************************************************
!BOP
! !IROUTINE: vmix_coeffs_ddiff
! !INTERFACE:

  subroutine vmix_coeffs_ddiff(Vmix_vars, Vmix_ddiff_params)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for the double diffusion mixing
!  parameterizatiion.
!\\
!\\
!
! !USES:
!  only those used by entire module.

! !INPUT PARAMETERS:
    type(vmix_ddiff_params_type), intent(in) :: Vmix_ddiff_params

! !INPUT/OUTPUT PARAMETERS:
    type(vmix_data_type), intent(inout) :: Vmix_vars
!EOP
!BOC

    Vmix_vars%visc_iface = Vmix_ddiff_params%deleteme
    Vmix_vars%diff_iface = Vmix_ddiff_params%deleteme

!EOC
  end subroutine vmix_coeffs_ddiff

end module vmix_ddiff
