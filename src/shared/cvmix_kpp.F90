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
   use cvmix_put_get, only :         cvmix_put
!EOP

   implicit none
   private
   save

!BOP

! !PUBLIC MEMBER FUNCTIONS:

   public :: cvmix_init_kpp
   public :: cvmix_coeffs_kpp
   public :: cvmix_put_kpp
   public :: cvmix_kpp_compute_OBL_depth ! MNL: I think we want this public

   interface cvmix_put_kpp
     module procedure cvmix_put_kpp_real
    end interface cvmix_put_kpp

! !PUBLIC TYPES:

  ! cvmix_kpp_params_type contains the necessary parameters for KPP mixing
  type, public :: cvmix_kpp_params_type
      private
      real(cvmix_r8) :: Ri_crit
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

    call cvmix_put_kpp(CVmix_kpp_params, 'Ri_crit', 0.3_cvmix_r8)
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

    call cvmix_kpp_compute_OBL_depth(CVmix_vars)
    CVmix_vars%visc_iface = CVmix_kpp_params%Ri_crit
    CVmix_vars%diff_iface = CVmix_kpp_params%Ri_crit

!EOC
  end subroutine cvmix_coeffs_kpp

!BOP

! !IROUTINE: cvmix_put_kpp_real
! !INTERFACE:

  subroutine cvmix_put_kpp_real(CVmix_kpp_params, varname, val)

! !DESCRIPTION:
!  Write a real value into a cvmix\_kpp\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    real(cvmix_r8),   intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_kpp_params_type), intent(inout) :: CVmix_kpp_params
!EOP
!BOC

    select case (trim(varname))
      case ('Ri_crit')
        CVmix_kpp_params%Ri_crit = val
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end subroutine cvmix_put_kpp_real

!BOP

! !IROUTINE: cvmix_kpp_compute_OBL_depth
! !INTERFACE:

  subroutine cvmix_kpp_compute_OBL_depth(CVMix_vars)

! !DESCRIPTION:
!  Computes the depth of the ocean boundary layer (OBL) for a given column
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:

! !OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars

!EOP
!BOC

    ! Local variables
    real(cvmix_r8) :: lcl_obl_depth

    lcl_obl_depth = real(-70, cvmix_r8)
    call cvmix_put(CVmix_vars, 'OBL_depth', lcl_obl_depth)

!EOC

  end subroutine cvmix_kpp_compute_OBL_depth

end module cvmix_kpp
