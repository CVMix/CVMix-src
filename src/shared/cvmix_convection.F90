module cvmix_convection

!BOP
!\newpage
! !MODULE: cvmix_convection
!
! !DESCRIPTION:
!  This module contains routines to initialize the derived types needed for
!  specifying mixing coefficients to parameterize vertical convective mixing,
!  and to set the viscosity and diffusivity in gravitationally unstable
!  portions of the water column.
!\\
!\\

! !USES:
  use cvmix_kinds_and_types, only : cvmix_r8,                                 &
                                    cvmix_zero,                               &
                                    cvmix_one,                                &
                                    cvmix_data_type,                          &
                                    CVMIX_OVERWRITE_OLD_VAL,                  &
                                    CVMIX_SUM_OLD_AND_NEW_VALS,               &
                                    CVMIX_MAX_OLD_AND_NEW_VALS

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

  interface cvmix_coeffs_conv
    module procedure cvmix_coeffs_conv_low
    module procedure cvmix_coeffs_conv_wrap
  end interface cvmix_coeffs_conv

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
      ! Threshold for squared buoyancy frequency needed to trigger
      ! Brunt-Vaisala parameterization
      real(cvmix_r8) :: BVsqr_convect ! units: s^-2
      ! Flag for what to do with old values of CVmix_vars%[MTS]diff
      integer :: handle_old_vals
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
    type (cvmix_conv_params_type), optional, intent(inout) ::                 &
                                       CVmix_conv_params_user

! !INPUT PARAMETERS:
    real(cvmix_r8), intent(in) :: &
      convect_diff,      &! diffusivity to parameterize convection
      convect_visc        ! viscosity to parameterize convection
    logical,        intent(in), optional :: lBruntVaisala ! True => B-V mixing
    real(cvmix_r8), intent(in), optional :: BVsqr_convect ! B-V parameter

!EOP
!BOC

    ! Set convect_diff and convect_visc in conv_params_type
    call cvmix_put_conv("convect_diff", convect_diff, CVmix_conv_params_user)
    call cvmix_put_conv("convect_visc", convect_visc, CVmix_conv_params_user)

    if (present(lBruntVaisala)) then
      call cvmix_put_conv("lBruntVaisala", lBruntVaisala,                     &
                          CVmix_conv_params_user)
    else
      call cvmix_put_conv("lBruntVaisala", .false., CVmix_conv_params_user)
    end if

    if (present(BVsqr_convect)) then
      call cvmix_put_conv("BVsqr_convect", BVsqr_convect,                     &
                          CVmix_conv_params_user)
    else
      call cvmix_put_conv("BVsqr_convect", cvmix_zero, CVmix_conv_params_user)
    end if

!EOC

  end subroutine cvmix_init_conv

!BOP

! !IROUTINE: cvmix_coeffs_conv_wrap
! !INTERFACE:

  subroutine cvmix_coeffs_conv_wrap(CVmix_vars, CVmix_conv_params_user)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for convective mixing.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:

    type (cvmix_conv_params_type), optional, target, intent(in)  ::           &
                                             CVmix_conv_params_user

! !INPUT/OUTPUT PARAMETERS:
    type (cvmix_data_type), intent(inout) :: CVmix_vars

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    real(cvmix_r8), dimension(:), allocatable :: new_Mdiff, new_Tdiff
    integer :: nlev, kw
    type (cvmix_conv_params_type), pointer :: CVmix_conv_params_in

    if (present(CVmix_conv_params_user)) then
      CVmix_conv_params_in => CVmix_conv_params_user
    else
      CVmix_conv_params_in => CVmix_conv_params_saved
    end if
    nlev = CVmix_vars%nlev
    allocate(new_Mdiff(nlev+1), new_Tdiff(nlev+1))

    call cvmix_coeffs_conv(new_Mdiff, new_Tdiff,                              &
                           CVmix_vars%SqrBuoyancyFreq_iface,                  &
                           CVmix_vars%WaterDensity_cntr,                      &
                           CVmix_vars%AdiabWaterDensity_cntr,                 &
                           CVmix_conv_params_user)

    select case (CVmix_conv_params_in%handle_old_vals)
      case (CVMIX_SUM_OLD_AND_NEW_VALS)
        call cvmix_put(CVmix_vars,"Mdiff", new_Mdiff + CVmix_vars%Mdiff_iface)
        call cvmix_put(CVmix_vars,"Tdiff", new_Tdiff + CVmix_vars%Tdiff_iface)
      case (CVMIX_MAX_OLD_AND_NEW_VALS)
        do kw=1,nlev+1
          CVmix_vars%Mdiff_iface(kw) = max(CVmix_vars%Mdiff_iface(kw),        &
                                           new_Mdiff(kw))
          CVmix_vars%Tdiff_iface(kw) = max(CVmix_vars%Tdiff_iface(kw),        &
                                           new_Tdiff(kw))
        end do
      case (CVMIX_OVERWRITE_OLD_VAL)
        call cvmix_put(CVmix_vars,"Mdiff", new_Mdiff)
        call cvmix_put(CVmix_vars,"Tdiff", new_Tdiff)
      case DEFAULT
        print*, "ERROR: do not know how to handle old values!"
        stop 1
    end select

    deallocate(new_Mdiff, new_Tdiff)

!EOC

  end subroutine cvmix_coeffs_conv_wrap
!BOP

! !IROUTINE: cvmix_coeffs_conv_low
! !INTERFACE:

  subroutine cvmix_coeffs_conv_low(Mdiff_out, Tdiff_out, Nsqr, dens,          &
                                   dens_lwr, CVmix_conv_params_user)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for convective mixing.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:

    ! nlev+1
    real(cvmix_r8), dimension(:), intent(in) :: Nsqr
    ! nlev
    real(cvmix_r8), dimension(:), intent(in) :: dens, dens_lwr
    type (cvmix_conv_params_type), optional, target, intent(in)  ::           &
                                             CVmix_conv_params_user

! !INPUT/OUTPUT PARAMETERS:
    ! nlev+1
    real(cvmix_r8), dimension(:), intent(inout) :: Mdiff_out, Tdiff_out

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    real(cvmix_r8) :: vvconv, wgt
    integer        :: nlev, kw
    type (cvmix_conv_params_type), pointer :: CVmix_conv_params_in

    if (present(CVmix_conv_params_user)) then
      CVmix_conv_params_in => CVmix_conv_params_user
    else
      CVmix_conv_params_in => CVmix_conv_params_saved
    end if
    nlev = size(dens)

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
        do kw=1, nlev
          wgt = cvmix_zero
          if (Nsqr(kw).le.0) then
            if (Nsqr(kw).gt.CVmix_conv_params_in%BVsqr_convect) then
              wgt = cvmix_one - Nsqr(kw) / CVmix_conv_params_in%BVsqr_convect
              wgt = (cvmix_one - wgt**2)**3
            else
              wgt = cvmix_one
            end if
          end if
          Mdiff_out(kw) = wgt*cvmix_get_conv_real('convect_visc',             &
                                                  CVmix_conv_params_in)
          Tdiff_out(kw) = wgt*cvmix_get_conv_real('convect_diff',             &
                                                  CVmix_conv_params_in)
        end do
      else ! BVsqr_convect >= 0 => step function
        do kw=1,nlev-1
          if (Nsqr(kw).le.0) then
            Mdiff_out(kw) = cvmix_get_conv_real('convect_visc',               &
                                                CVmix_conv_params_in)
            Tdiff_out(kw) = cvmix_get_conv_real('convect_diff',               &
                                                CVmix_conv_params_in)
          else
            Mdiff_out(kw) = cvmix_zero
            Tdiff_out(kw) = cvmix_zero
          end if
        end do
        Mdiff_out(nlev) = cvmix_zero
        Tdiff_out(nlev) = cvmix_zero
      end if
      Mdiff_out(nlev+1) = cvmix_zero
      Tdiff_out(nlev+1) = cvmix_zero
    else
      ! Default convection mixing based on density
      do kw=1,nlev-1
        if (CVmix_conv_params_in%convect_visc.ne.cvmix_zero) then
           vvconv = cvmix_get_conv_real('convect_visc', CVmix_conv_params_in)
        else
          ! convection only affects tracers
          vvconv = Mdiff_out(kw)
        end if

!        if (CVmix_vars%WaterDensity_cntr(kw).gt.                              &
!            CVmix_vars%AdiabWaterDensity_cntr(kw)) then
        if (dens(kw).gt.dens_lwr(kw)) then
          Mdiff_out(kw+1) = vvconv
          Tdiff_out(kw+1) = cvmix_get_conv_real('convect_diff',               &
                                                CVmix_conv_params_in)
        end if
      end do
    end if

!EOC

  end subroutine cvmix_coeffs_conv_low

!BOP

! !IROUTINE: cvmix_put_conv_real
! !INTERFACE:

  subroutine cvmix_put_conv_real(varname, val, CVmix_conv_params_user)

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
    type (cvmix_conv_params_type), optional, target, intent(inout) ::         &
                                               CVmix_conv_params_user

!EOP
!BOC

    type (cvmix_conv_params_type), pointer :: CVmix_conv_params_out

    if (present(CVmix_conv_params_user)) then
      CVmix_conv_params_out => CVmix_conv_params_user
    else
      CVmix_conv_params_out => CVmix_conv_params_saved
    end if

    select case (trim(varname))
      case ('convect_diff')
        CVmix_conv_params_out%convect_diff = val
      case ('convect_visc')
        CVmix_conv_params_out%convect_visc = val
      case ('BVsqr_convect')
        CVmix_conv_params_out%BVsqr_convect = val
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end subroutine cvmix_put_conv_real

!BOP

! !IROUTINE: cvmix_put_conv_logical
! !INTERFACE:

  subroutine cvmix_put_conv_logical(varname, val, CVmix_conv_params_user)

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
    type (cvmix_conv_params_type), optional, target, intent(inout) ::         &
                                               CVmix_conv_params_user

!EOP
!BOC

    type (cvmix_conv_params_type), pointer :: CVmix_conv_params_out

    if (present(CVmix_conv_params_user)) then
      CVmix_conv_params_out => CVmix_conv_params_user
    else
      CVmix_conv_params_out => CVmix_conv_params_saved
    end if

    select case (trim(varname))
      case ('lBruntVaisala')
        CVmix_conv_params_out%lBruntVaisala = val
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

    cvmix_get_conv_real = cvmix_zero
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

