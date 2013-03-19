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
                                     cvmix_data_type,        &
                                     cvmix_conv_params_type
   use cvmix_put_get, only : cvmix_put
!EOP

  implicit none
  private
  save

!BOP

! !PUBLIC MEMBER FUNCTIONS:

   public :: cvmix_init_conv
   public :: cvmix_coeffs_conv
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
    call cvmix_put(CVmix_conv_params, "convect_diff", convect_diff)
    call cvmix_put(CVmix_conv_params, "convect_visc", convect_visc)

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

end module cvmix_convection

