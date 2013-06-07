module cvmix_convection

!BOP
!\newpage
! !MODULE: cvmix_convection
!
! !DESCRIPTION:
!  This module contains routines to initialize the derived types needed for
!  specifying mixing coefficients to parameterize vertical convective mixing,
!  and to set the viscosity and diffusivity in gravitationally unstable portions
!  of the water column.
!\\
!\\

! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !USES:
   use cvmix_kinds_and_types, only : cvmix_r8,               &
                                     cvmix_data_type
!EOP

  implicit none
  private
  save

!BOP

! !PUBLIC MEMBER FUNCTIONS:

   public :: cvmix_init_conv
   public :: cvmix_coeffs_conv
   public :: cvmix_put_conv
   public :: cvmix_get_conv_real

   interface cvmix_put_conv
     module procedure cvmix_put_conv_real
   end interface cvmix_put_conv

  ! cvmix_conv_params_type contains the necessary parameters for convective
  ! mixing.
  type, public :: cvmix_conv_params_type
    private
    real(cvmix_r8) :: convect_diff
    real(cvmix_r8) :: convect_visc
  end type cvmix_conv_params_type
!EOP

contains

!BOP

! !IROUTINE: cvmix_init_conv
! !INTERFACE:

  subroutine cvmix_init_conv(CVmix_conv_params, convect_diff, convect_visc)

! !DESCRIPTION:
!  Initialization routine for specifying convective mixing coefficients.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !OUTPUT PARAMETERS:
    type (cvmix_conv_params_type), intent(out) :: CVmix_conv_params

! !INPUT PARAMETERS:
   real(cvmix_r8), intent(in) :: &
      convect_diff,      &! diffusivity to parameterize convection
      convect_visc        ! viscosity to parameterize convection
!EOP
!BOC

    ! Set convect_diff and convect_visc in conv_params_type
    call cvmix_put_conv(CVmix_conv_params, "convect_diff", convect_diff)
    call cvmix_put_conv(CVmix_conv_params, "convect_visc", convect_visc)

!EOC
  end subroutine cvmix_init_conv


!BOP

! !IROUTINE: cvmix_coeffs_conv
! !INTERFACE:

  subroutine cvmix_coeffs_conv(CVmix_vars, CVmix_conv_params)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for convective mixing.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:

    type (cvmix_conv_params_type), intent(in)  :: CVmix_conv_params

! !INPUT/OUTPUT PARAMETERS:
    type (cvmix_data_type), intent(inout) :: CVmix_vars
!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    real(cvmix_r8) :: vvconv
    integer        :: kw  ! vertical level index 

!-----------------------------------------------------------------------
!
!  enhance the vertical mixing coefficients if gravitationally unstable
!
!-----------------------------------------------------------------------

    do kw=1,CVmix_vars%nlev-1
      if (CVmix_conv_params%convect_visc.ne.0_cvmix_r8) then
         vvconv = CVmix_conv_params%convect_visc
      else
        ! convection only affects tracers
        vvconv = CVmix_vars%visc_iface(kw)
      end if

      if (CVmix_vars%dens(kw).gt.CVmix_vars%dens_lwr(kw)) then
        CVmix_vars%diff_iface(kw+1,1) = CVmix_conv_params%convect_diff
        CVmix_vars%visc_iface(kw+1)   = vvconv
      end if
    end do

!EOC

  end subroutine cvmix_coeffs_conv

!BOP

! !IROUTINE: cvmix_put_conv_real
! !INTERFACE:

  subroutine cvmix_put_conv_real(CVmix_conv_params, varname, val)

! !DESCRIPTION:
!  Write a real value into a cvmix\_conv\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    real(cvmix_r8),   intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_conv_params_type), intent(inout) :: CVmix_conv_params
!EOP
!BOC

    select case (trim(varname))
      case ('convect_diff')
        CVmix_conv_params%convect_diff = val
      case ('convect_visc')
        CVmix_conv_params%convect_visc = val
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end subroutine cvmix_put_conv_real

!BOP

! !IROUTINE: cvmix_get_conv_real
! !INTERFACE:

  function cvmix_get_conv_real(CVmix_conv_params, varname)

! !DESCRIPTION:
!  Read the real value of a cvmix\_conv\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    type(cvmix_conv_params_type), intent(in) :: CVmix_conv_params
    character(len=*),             intent(in) :: varname

! !OUTPUT PARAMETERS:
    real(cvmix_r8) :: cvmix_get_conv_real
!EOP
!BOC

    select case (trim(varname))
      case ('convect_diff')
        cvmix_get_conv_real = CVmix_conv_params%convect_diff
      case ('convect_visc')
        cvmix_get_conv_real = CVmix_conv_params%convect_visc
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end function cvmix_get_conv_real

end module cvmix_convection

