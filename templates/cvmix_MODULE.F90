module cvmix_MODULE

!BOP
!\newpage
! !MODULE: cvmix_MODULE
!
! !DESCRIPTION:
!  This module contains routines to initialize the derived types needed for
!  MODULE_MIX_TYPE mixing and to set the diffusivity coefficient
!  accordingly.
!\\
!\\
!
! !REVISION HISTORY:
!  SVN:$Id: cvmix_ddiff.F90 258 2013-10-21 15:17:26Z mike.levy.work@gmail.com $
!  SVN:$URL: https://cvmix.googlecode.com/svn/trunk/src/shared/cvmix_ddiff.F90 $

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                                 &
                                    cvmix_data_type
!EOP

  implicit none
  private
  save

!BOP

! !PUBLIC MEMBER FUNCTIONS:

  public :: cvmix_init_MODULE
  public :: cvmix_coeffs_MODULE
  public :: cvmix_put_MODULE
  public :: cvmix_get_MODULE_real

  interface cvmix_put_MODULE
    module procedure cvmix_put_MODULE_real
  end interface cvmix_put_MODULE

! !PUBLIC TYPES:

  ! cvmix_MODULE_params_type contains the necessary parameters for
  ! MODULE_MIX_TYPE mixing
  type, public :: cvmix_MODULE_params_type
    private
    real(cvmix_r8) :: VARNAME
  end type cvmix_MODULE_params_type

!EOP

  type(cvmix_MODULE_params_type), target :: CVmix_MODULE_params_saved

 contains

!BOP

! !IROUTINE: cvmix_init_MODULE
! !INTERFACE:

  subroutine cvmix_init_MODULE(VARNAME, CVmix_MODULE_params_user)

! !DESCRIPTION:
!  Initialization routine for MODULE_MIX_TYPE mixing. This mixing technique...
!  (description of technique, use latex formatting)
!\\
! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
! These should all be optional, with values set to default if variable is not
! present
  real(cvmix_r8), optional :: VARNAME

! !OUTPUT PARAMETERS:
    type(cvmix_MODULE_params_type), optional, target, intent(inout) ::        &
                                              CVmix_MODULE_params_user
!EOP
!BOC

    type(cvmix_MODULE_params_type), pointer :: CVmix_MODULE_params_out

    if (present(CVmix_MODULE_params_user)) then
      CVmix_MODULE_params_out => CVmix_MODULE_params_user
    else
      CVmix_MODULE_params_out => CVmix_MODULE_params_saved
    end if

    ! For each variable, check if present and either set to initialized
    ! value or to default value
    if (present(VARNAME)) then
      ! If VARNAME is present, set CVmix_MODULE_params_out%VARNAME = VARNAME
      call cvmix_put_MODULE(CVmix_MODULE_params_out, "VARNAME", VARNAME)
    else
      ! Else set to default value (0 in this case)
      call cvmix_put_MODULE(CVmix_MODULE_params_out, "VARNAME", 0.0_cvmix_r8)
    end if

!EOC

  end subroutine cvmix_init_MODULE

!MNL NOTE: Still need to break ddiff into low- / wrapped- interfaces.

!BOP

! !IROUTINE: cvmix_put_MODULE_real
! !INTERFACE:

  subroutine cvmix_put_MODULE_real(CVmix_MODULE_params, varname, val)

! !DESCRIPTION:
!  Write a real value into a cvmix\_MODULE\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    real(cvmix_r8),   intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_MODULE_params_type), intent(inout) :: CVmix_MODULE_params
!EOP
!BOC

    select case (trim(varname))
      case ('VARNAME')
        CVmix_MODULE_params%VARNAME = val
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end subroutine cvmix_put_MODULE_real

!BOP

! !IROUTINE: cvmix_get_VARNAME_real
! !INTERFACE:

  function cvmix_get_VARNAME_real(varname, CVmix_MODULE_params_user)

! !DESCRIPTION:
!  Return the real value of a cvmix\_MODULE\_params\_type variable.
!  NOTE: This function is not efficient and is only for infrequent
!  queries of MODULE parameters, such as at initialization.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*),                                intent(in)    :: varname
    type(cvmix_MODULE_params_type), optional, target, intent(inout) ::        &
                                              CVmix_MODULE_params_user

! !OUTPUT PARAMETERS:
    real(cvmix_r8) :: cvmix_get_MODULE_real
!EOP
!BOC

    type(cvmix_MODULE_params_type), pointer :: CVmix_MODULE_params_get

    if (present(CVmix_MODULE_params_user)) then
      CVmix_MODULE_params_get => CVmix_MODULE_params_user
    else
      CVmix_MODULE_params_get => CVmix_MODULE_params_saved
    end if

    cvmix_get_MODULE_real = 0.0_cvmix_r8
    select case (trim(varname))
      case ('VARNAME')
        cvmix_get_MODULE_real = CVmix_MODULE_params_get%VARNAME
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end function cvmix_get_MODULE_real

end module cvmix_MODULE
