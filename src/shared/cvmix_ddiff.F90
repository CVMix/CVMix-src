!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module cvmix_ddiff

!BOP
!\newpage
! !MODULE: cvmix_ddiff
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

   use cvmix_kinds_and_types, only : cvmix_r8,                &
                                     cvmix_data_type,         &
                                     cvmix_ddiff_params_type
!EOP

   implicit none
   private
   save

!BOP

! !PUBLIC MEMBER FUNCTIONS:

   public :: cvmix_init_ddiff
   public :: cvmix_coeffs_ddiff
!EOP

 contains

!BOP

! !IROUTINE: cvmix_init_ddiff
! !INTERFACE:

  subroutine cvmix_init_ddiff(CVmix_ddiff_params)

! !DESCRIPTION:
!  Initialization routine for double diffusion mixing.
!
! !USES:
!  Only those used by entire module.

! !OUTPUT PARAMETERS:
    type(cvmix_ddiff_params_type), intent(inout) :: CVmix_ddiff_params
!EOP
!BOC

    CVmix_ddiff_params%deleteme = 0.0_cvmix_r8
!EOC

  end subroutine cvmix_init_ddiff

!***********************************************************************
!BOP
! !IROUTINE: cvmix_coeffs_ddiff
! !INTERFACE:

  subroutine cvmix_coeffs_ddiff(CVmix_vars, CVmix_ddiff_params)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for the double diffusion mixing
!  parameterizatiion.
!\\
!\\
!
! !USES:
!  only those used by entire module.

! !INPUT PARAMETERS:
    type(cvmix_ddiff_params_type), intent(in) :: CVmix_ddiff_params

! !INPUT/OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars
!EOP
!BOC

    CVmix_vars%visc_iface = CVmix_ddiff_params%deleteme
    CVmix_vars%diff_iface = CVmix_ddiff_params%deleteme

!EOC
  end subroutine cvmix_coeffs_ddiff

end module cvmix_ddiff
