module vmix_convection

!BOP
! !MODULE: vmix_convection
!
! !DESCRIPTION:
!  This module contains routines to initialize the derived types needed for
!  convective mixing and to set the viscosity and diffusivity in unstable
!  columns.
!\\
!\\

! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !USES:
   use vmix_kinds_and_types, only : vmix_r8,               &
                                    vmix_data_type,        &
                                    vmix_conv_params_type

!EOP

  implicit none
  private
  save

!BOP

! !PUBLIC MEMBER FUNCTIONS:

   public :: vmix_init_conv
   public :: vmix_coeffs_conv

!EOP

contains

!BOP

! !IROUTINE: vmix_init_conv
! !INTERFACE:

  subroutine vmix_init_conv(Vmix_conv_params, convect_diff, convect_visc)

! !DESCRIPTION:
!  Initialization routine for convective mixing.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !OUTPUT PARAMETERS:
    type (vmix_conv_params_type), intent(out) :: Vmix_conv_params

! !INPUT PARAMETERS:
   real(vmix_r8), intent(in) :: &
      convect_diff,      &! diffusivity to mimic convection
      convect_visc        ! viscosity to mimic convection

!EOP
!BOC

    ! Set convect_diff and convect_visc in conv_params_type
    Vmix_conv_params%convect_diff = convect_diff
    Vmix_conv_params%convect_visc = convect_visc

!EOC
  end subroutine vmix_init_conv


!BOP

! !IROUTINE: vmix_coeffs_conv
! !INTERFACE:

  subroutine vmix_coeffs_conv(Vmix_vars, Vmix_conv_params)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for convective mixing.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:

    type (vmix_conv_params_type), intent(in)  :: Vmix_conv_params

! !INPUT/OUTPUT PARAMETERS:
    type (vmix_data_type), intent(inout) :: Vmix_vars


!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer :: kw  ! vertical level index 
       
    real(vmix_r8) :: vvconv

!-----------------------------------------------------------------------
!
!  enhance the vertical mixing coefficients if gravitationally unstable
!
!-----------------------------------------------------------------------

    do kw=1,Vmix_vars%nlev-1
      if (Vmix_conv_params%convect_visc.ne.0_vmix_r8) then
         vvconv = Vmix_conv_params%convect_visc
      else
        ! convection only affects tracers
        vvconv = Vmix_vars%visc_iface(kw)
      end if

      if (Vmix_vars%dens(kw).gt.Vmix_vars%dens_lwr(kw)) then
        Vmix_vars%diff_iface(kw+1,1) = Vmix_conv_params%convect_diff
        Vmix_vars%visc_iface(kw+1)   = vvconv
      end if
    end do

!EOC

  end subroutine vmix_coeffs_conv

end module vmix_convection

