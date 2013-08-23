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
     module procedure cvmix_put_conv_logical
   end interface cvmix_put_conv

! !PUBLIC TYPES:

  ! cvmix_conv_params_type contains the necessary parameters for convective
  ! mixing.
  type, public :: cvmix_conv_params_type
    private
    ! Convective diff
    ! diffusivity coefficient used in convective regime
    real(cvmix_r8) :: convect_diff ! units: m^2/s
    ! viscosity coefficient used in convective regime
    real(cvmix_r8) :: convect_visc ! units: m^2/s
    logical        :: lBruntVaisala
    ! Threshold for squared buoyancy frequency needed to trigger Brunt-Vaisala
    ! parameterization
    real(cvmix_r8) :: BVsqr_convect ! units: s^-2
  end type cvmix_conv_params_type

!EOP

  type(cvmix_conv_params_type), target :: CVmix_conv_params_saved

contains

!BOP

! !IROUTINE: cvmix_init_conv
! !INTERFACE:

  subroutine cvmix_init_conv(convect_diff, convect_visc, lBruntVaisala,       &
                             BVsqr_convect, CVmix_conv_params_user)

! !DESCRIPTION:
!  Initialization routine for specifying convective mixing coefficients.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !OUTPUT PARAMETERS:
    type (cvmix_conv_params_type), optional, target, intent(out) ::           &
                                                        CVmix_conv_params_user

! !INPUT PARAMETERS:
    real(cvmix_r8), intent(in) :: &
      convect_diff,      &! diffusivity to parameterize convection
      convect_visc        ! viscosity to parameterize convection
    logical,        intent(in), optional :: lBruntVaisala ! True => B-V mixing
    real(cvmix_r8), intent(in), optional :: BVsqr_convect ! B-V parameter
!EOP
!BOC

    type (cvmix_conv_params_type), pointer :: CVmix_conv_params_out

    if (present(CVmix_conv_params_user)) then
      CVmix_conv_params_out => CVmix_conv_params_user
    else
      CVmix_conv_params_out => CVmix_conv_params_saved
    end if

    ! Set convect_diff and convect_visc in conv_params_type
    call cvmix_put_conv(CVmix_conv_params_out, "convect_diff", convect_diff)
    call cvmix_put_conv(CVmix_conv_params_out, "convect_visc", convect_visc)

    if (present(lBruntVaisala)) then
      call cvmix_put_conv(CVmix_conv_params_out, "lBruntVaisala",             &
                          lBruntVaisala)
    else
      call cvmix_put_conv(CVmix_conv_params_out, "lBruntVaisala", .false.)
    end if

    if (present(BVsqr_convect)) then
      call cvmix_put_conv(CVmix_conv_params_out, "BVsqr_convect",             &
                          BVsqr_convect)
    else
      call cvmix_put_conv(CVmix_conv_params_out, "BVsqr_convect", 0.0_cvmix_r8)
    end if

!EOC

  end subroutine cvmix_init_conv

!BOP

! !IROUTINE: cvmix_coeffs_conv
! !INTERFACE:

  subroutine cvmix_coeffs_conv(CVmix_vars, CVmix_conv_params_user)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for convective mixing.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:

    type (cvmix_conv_params_type), optional, target, intent(in)  :: CVmix_conv_params_user

! !INPUT/OUTPUT PARAMETERS:
    type (cvmix_data_type), intent(inout) :: CVmix_vars
!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    real(cvmix_r8) :: vvconv, wgt
    integer        :: kw  ! vertical level index 

    type (cvmix_conv_params_type), pointer :: CVmix_conv_params_in

    if (present(CVmix_conv_params_user)) then
      CVmix_conv_params_in => CVmix_conv_params_user
    else
      CVmix_conv_params_in => CVmix_conv_params_saved
    end if

!-----------------------------------------------------------------------
!
!  enhance the vertical mixing coefficients if gravitationally unstable
!
!-----------------------------------------------------------------------
    if (CVmix_conv_params_in%lBruntVaisala) then
      ! Brunt-Vaisala mixing based on buoyancy
      ! Based on parameter BVsqr_convect
      ! diffusivity = convect_diff * wgt
      ! viscosity   = convect_visc * wgt

      ! For BVsqr_convect < 0:
      ! wgt = 0 for N^2 > 0
      ! wgt = 1 for N^2 < BVsqr_convect
      ! wgt = [1 - (1-N^2/BVsqr_convect)^2]^3 otherwise

      ! If BVsqr_convect >= 0:
      ! wgt = 0 for N^2 > 0
      ! wgt = 1 for N^2 <= 0

      ! Compute wgt
      if (CVmix_conv_params_in%BVsqr_convect.lt.0) then
        do kw=1,CVmix_vars%nlev
          wgt = 0.0_cvmix_r8
          if (CVmix_vars%buoy_iface(kw).le.0) then
            if (CVmix_vars%buoy_iface(kw).gt.                                 &
                CVmix_conv_params_in%BVsqr_convect) then
              wgt = 1.0_cvmix_r8 - CVmix_vars%buoy_iface(kw) /                &
                    CVmix_conv_params_in%BVsqr_convect
              wgt = (1.0_cvmix_r8 - wgt**2)**3
            else
              wgt = 1.0_cvmix_r8
            end if
          end if
          CVmix_vars%visc_iface(kw) = wgt*cvmix_get_conv_real('convect_visc', &
                                          CVmix_conv_params_in)
          CVmix_vars%diff_iface(kw,:)=wgt*cvmix_get_conv_real('convect_diff', &
                                          CVmix_conv_params_in)
        end do
      else ! BVsqr_convect >= 0 => step function
        do kw=1,CVmix_vars%nlev-1
          if (CVmix_vars%buoy_iface(kw).le.0) then
            CVmix_vars%visc_iface(kw)   = cvmix_get_conv_real('convect_visc', &
                                          CVmix_conv_params_in)
            CVmix_vars%diff_iface(kw,:) = cvmix_get_conv_real('convect_diff', &
                                          CVmix_conv_params_in)
          else
            CVmix_vars%visc_iface(kw)   = 0.0_cvmix_r8
            CVmix_vars%diff_iface(kw,:) = 0.0_cvmix_r8
          end if
        end do
      end if
      CVmix_vars%visc_iface(CVmix_vars%nlev+1)   = 0.0_cvmix_r8
      CVmix_vars%diff_iface(CVmix_vars%nlev+1,:) = 0.0_cvmix_r8
    else
      ! Default convection mixing based on density
      do kw=1,CVmix_vars%nlev-1
        if (CVmix_conv_params_in%convect_visc.ne.0_cvmix_r8) then
           vvconv = cvmix_get_conv_real('convect_visc', CVmix_conv_params_in)
        else
          ! convection only affects tracers
          vvconv = CVmix_vars%visc_iface(kw)
        end if

        if (CVmix_vars%dens(kw).gt.CVmix_vars%dens_lwr(kw)) then
          CVmix_vars%diff_iface(kw+1,1) = cvmix_get_conv_real('convect_diff', &
                                          CVmix_conv_params_in)
          CVmix_vars%visc_iface(kw+1)   = vvconv
        end if
      end do
    end if

!EOC

  end subroutine cvmix_coeffs_conv

!BOP

! !IROUTINE: cvmix_put_conv_real
! !INTERFACE:

  subroutine cvmix_put_conv_real(CVmix_conv_params_put, varname, val)

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
    type(cvmix_conv_params_type), intent(inout) :: CVmix_conv_params_put
!EOP
!BOC

    select case (trim(varname))
      case ('convect_diff')
        CVmix_conv_params_put%convect_diff = val
      case ('convect_visc')
        CVmix_conv_params_put%convect_visc = val
      case ('BVsqr_convect')
        CVmix_conv_params_put%BVsqr_convect = val
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end subroutine cvmix_put_conv_real

!BOP

! !IROUTINE: cvmix_put_conv_logical
! !INTERFACE:

  subroutine cvmix_put_conv_logical(CVmix_conv_params_put, varname, val)

! !DESCRIPTION:
!  Write a Boolean value into a cvmix\_conv\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    logical,          intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_conv_params_type), intent(inout) :: CVmix_conv_params_put
!EOP
!BOC

    select case (trim(varname))
      case ('lBruntVaisala')
        CVmix_conv_params_put%lBruntVaisala = val
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
    end select

!EOC

  end subroutine cvmix_put_conv_logical

!BOP

! !IROUTINE: cvmix_get_conv_real
! !INTERFACE:

  function cvmix_get_conv_real(varname, CVmix_conv_params_user)

! !DESCRIPTION:
!  Read the real value of a cvmix\_conv\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*),             intent(in) :: varname
    type(cvmix_conv_params_type), optional, target, intent(in) ::             &
                                           CVmix_conv_params_user

! !OUTPUT PARAMETERS:
    real(cvmix_r8) :: cvmix_get_conv_real
!EOP
!BOC

    type(cvmix_conv_params_type), pointer :: CVmix_conv_params_get

    if (present(CVmix_conv_params_user)) then
      CVmix_conv_params_get => CVmix_conv_params_user
    else
      CVmix_conv_params_get => CVmix_conv_params_saved
    end if

    cvmix_get_conv_real = 0.0_cvmix_r8
    select case (trim(varname))
      case ('convect_diff')
        cvmix_get_conv_real = CVmix_conv_params_get%convect_diff
      case ('convect_visc')
        cvmix_get_conv_real = CVmix_conv_params_get%convect_visc
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end function cvmix_get_conv_real

end module cvmix_convection

