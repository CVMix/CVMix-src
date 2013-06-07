!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module cvmix_kpp

!BOP
!\newpage
! !MODULE: cvmix_kpp
!
! !DESCRIPTION:
!  This module contains routines to initialize the derived types needed for
!  KPP mixing and to set the viscosity and diffusivity coefficients
!  accordingly.
!\\
!\\
!
! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !USES:

   use cvmix_kinds_and_types, only : cvmix_r8,                &
                                     cvmix_data_type
!EOP

   implicit none
   private
   save

!BOP

! !PUBLIC MEMBER FUNCTIONS:

   public :: cvmix_init_kpp
   public :: cvmix_coeffs_kpp

! !PUBLIC TYPES:

  ! cvmix_kpp_params_type contains the necessary parameters for KPP mixing
  type, public :: cvmix_kpp_params_type
      private
      real(cvmix_r8) :: deleteme
  end type cvmix_kpp_params_type

!EOP

 contains

!BOP

! !IROUTINE: cvmix_init_kpp
! !INTERFACE:

  subroutine cvmix_init_kpp(CVmix_kpp_params)

! !DESCRIPTION:
!  Initialization routine for KPP mixing.
!
! !USES:
!  Only those used by entire module.

! !OUTPUT PARAMETERS:
    type(cvmix_kpp_params_type), intent(inout) :: CVmix_kpp_params

!EOP
!BOC

    CVmix_kpp_params%deleteme = 0.0_cvmix_r8
!EOC

  end subroutine cvmix_init_kpp

!***********************************************************************
!BOP
! !IROUTINE: cvmix_coeffs_kpp
! !INTERFACE:

  subroutine cvmix_coeffs_kpp(CVmix_vars, CVmix_kpp_params)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for the double diffusion mixing
!  parameterizatiion.
!\\
!\\
!
! !USES:
!  only those used by entire module.

! !INPUT PARAMETERS:
    type(cvmix_kpp_params_type), intent(in) :: CVmix_kpp_params

! !INPUT/OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars

!EOP
!BOC

    CVmix_vars%visc_iface = CVmix_kpp_params%deleteme
    CVmix_vars%diff_iface = CVmix_kpp_params%deleteme

!EOC
  end subroutine cvmix_coeffs_kpp

end module cvmix_kpp
