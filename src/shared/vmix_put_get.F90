module vmix_put_get

!BOP
!\newpage
! !MODULE: vmix_put_get
!
! !DESCRIPTION:
!  This module contains routines to pack data into the vmix datatypes
!  (allocating memory as necessary) and then unpack the data out. If we switch
!  to pointers, the pack will just point at the right target and the unpack
!  will be un-necessary.
!\\
!\\

! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !USES:

   use vmix_kinds_and_types, only : vmix_r8,                  &
                                    vmix_strlen,              &
                                    vmix_data_type,           &
                                    vmix_global_params_type,  &
                                    vmix_bkgnd_params_type
!EOP

  implicit none
  private
  save

!BOP

! !PUBLIC MEMBER FUNCTIONS:
   public :: vmix_put

  interface vmix_put
    module procedure vmix_put_int
    module procedure vmix_put_real
    module procedure vmix_put_real_1D
    module procedure vmix_put_bkgnd_real_1D
    module procedure vmix_put_global_params_int
    module procedure vmix_put_global_params_real
  end interface vmix_put
!EOP
!BOC

contains

!BOP

! !IROUTINE: vmix_put_int
! !INTERFACE:

  subroutine vmix_put_int(Vmix_vars, varname, val, opts)

! !DESCRIPTION:
!  Write an integer value into a vmix\_data\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*),           intent(in) :: varname
    integer,                    intent(in) :: val
    character(len=*), optional, intent(in) :: opts

! !OUTPUT PARAMETERS:
    type (vmix_data_type), intent(inout) :: Vmix_vars
!EOP
!BOC

    ! Local variables
    integer :: nlev
  
    nlev = Vmix_vars%nlev

    if ((trim(varname).ne.'nlev').and.(nlev.eq.-1)) then
      print*, "ERROR: you must specify the number of levels before ", &
              "you can pack data into a vmix_data_type!"
      print*, "You tried to set ", trim(varname)
      stop
    end if
    
    select case (trim(varname))
      case ('nlev')
        Vmix_vars%nlev = val
        
      case ('diff')
      if (.not.associated(Vmix_vars%diff_iface)) then
        allocate(Vmix_vars%diff_iface(nlev+1,2))
      end if
      if (present(opts)) then
        if (trim(opts).eq.'col2') then
          Vmix_vars%diff_iface(:,2) = val
        else
          Vmix_vars%diff_iface(:,1) = val
        end if
      end if
      
      case ('visc')
      if (.not.associated(Vmix_vars%visc_iface)) then
        allocate(Vmix_vars%visc_iface(nlev+1))
      end if
      Vmix_vars%visc_iface(:) = val

      case default
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop
      
    end select
!EOC

  end subroutine vmix_put_int

!BOP

! !IROUTINE: vmix_put_real
! !INTERFACE:

  subroutine vmix_put_real(Vmix_vars, varname, val, opts)

! !DESCRIPTION:
!  Write a real value into a vmix\_data\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*),           intent(in) :: varname
    real(vmix_r8),              intent(in) :: val
    character(len=*), optional, intent(in) :: opts

! !OUTPUT PARAMETERS:
    type (vmix_data_type), intent(inout) :: Vmix_vars
!EOP
!BOC

    ! Local variables
    integer :: nlev
  
    nlev = Vmix_vars%nlev

    if (nlev.eq.-1) then
      print*, "ERROR: you must specify the number of levels before ", &
              "you can pack data into a vmix_data_type!"
      print*, "You tried to set ", trim(varname)
      stop
    end if
    
    select case (trim(varname))
      case ('diff')
      if (.not.associated(Vmix_vars%diff_iface)) then
        allocate(Vmix_vars%diff_iface(nlev+1,2))
      end if
      if (present(opts)) then
        if (trim(opts).eq.'col2') then
          Vmix_vars%diff_iface(:,2) = val
        else
          Vmix_vars%diff_iface(:,1) = val
        end if
      end if
      
      case ('visc')
      if (.not.associated(Vmix_vars%visc_iface)) then
        allocate(Vmix_vars%visc_iface(nlev+1))
      end if
      Vmix_vars%visc_iface(:) = val

      case default
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop
      
    end select
!EOC

  end subroutine vmix_put_real

!BOP

! !IROUTINE: vmix_put_real_1D
! !INTERFACE:

  subroutine vmix_put_real_1D(Vmix_vars, varname, val, opts)

! !DESCRIPTION:
!  Write an array of real values into a vmix\_data\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*),            intent(in) :: varname
    real(vmix_r8), dimension(:), intent(in) :: val
    character(len=*), optional,  intent(in) :: opts

! !OUTPUT PARAMETERS:
    type (vmix_data_type), intent(inout) :: Vmix_vars
!EOP
!BOC

    ! Local variables
    integer :: nlev
  
    nlev = Vmix_vars%nlev

    if (nlev.eq.-1) then
      print*, "ERROR: you must specify the number of levels before ", &
              "you can pack data into a vmix_data_type!"
      print*, "You tried to set ", trim(varname)
      stop
    end if
    
    select case (trim(varname))
      case ('diff')
      if (.not.associated(Vmix_vars%diff_iface)) then
        allocate(Vmix_vars%diff_iface(nlev+1,2))
      end if
      if (present(opts)) then
        if (trim(opts).eq.'col2') then
          Vmix_vars%diff_iface(:,2) = val
        else
          Vmix_vars%diff_iface(:,1) = val
        end if
      end if

      case ('visc')
      if (.not.associated(Vmix_vars%visc_iface)) then
        allocate(Vmix_vars%visc_iface(nlev+1))
      end if
      Vmix_vars%visc_iface(:) = val

      case ('dens')
      if (.not.associated(Vmix_vars%dens)) then
        allocate(Vmix_vars%dens(nlev))
      end if
      Vmix_vars%dens(:) = val

      case ('dens_lwr')
      if (.not.associated(Vmix_vars%dens_lwr)) then
        allocate(Vmix_vars%dens_lwr(nlev))
      end if
      Vmix_vars%dens_lwr(:) = val

      case ('zw')
      if (.not.associated(Vmix_vars%z_iface)) then
        allocate(Vmix_vars%z_iface(nlev+1))
      end if
      Vmix_vars%z_iface(:) = val

      case default
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop
      
    end select

!EOC

  end subroutine vmix_put_real_1D

!BOP

! !IROUTINE: vmix_put_bkgnd_real_1D
! !INTERFACE:

  subroutine vmix_put_bkgnd_real_1D(varname, Vmix_bkgnd_params, val, & 
                                    opts)

! !DESCRIPTION:
!  Write an array of real values into a vmix\_bkgnd\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=vmix_strlen),  intent(in)          :: varname
    real(vmix_r8), dimension(:), intent(in)          :: val
    character(len=vmix_strlen), optional, intent(in) :: opts

! !OUTPUT PARAMETERS:
    type (vmix_bkgnd_params_type), intent(out) :: Vmix_bkgnd_params
!EOP
!BOC

!EOC

  end subroutine vmix_put_bkgnd_real_1D

!BOP

! !IROUTINE: vmix_put_global_params_int
! !INTERFACE:

  subroutine vmix_put_global_params_int(Vmix_params, varname, val, opts)

! !DESCRIPTION:
!  Write an integer value into a vmix\_global\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*),           intent(in) :: varname
    integer,                    intent(in) :: val
    character(len=*), optional, intent(in) :: opts

! !OUTPUT PARAMETERS:
    type (vmix_global_params_type), intent(inout) :: Vmix_params
!EOP
!BOC

    select case (trim(varname))
      case ('max_nlev')
        Vmix_params%max_nlev = val
        
      case default
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop
      
    end select
!EOC

  end subroutine vmix_put_global_params_int

!BOP

! !IROUTINE: vmix_put_global_params_real
! !INTERFACE:

  subroutine vmix_put_global_params_real(Vmix_params, varname, val, opts)

! !DESCRIPTION:
!  Write a real value into a vmix\_global\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*),           intent(in) :: varname
    real(vmix_r8),              intent(in) :: val
    character(len=*), optional, intent(in) :: opts

! !OUTPUT PARAMETERS:
    type (vmix_global_params_type), intent(inout) :: Vmix_params
!EOP
!BOC

    select case (trim(varname))
      case ('prandtl')
        Vmix_params%prandtl = val
        
      case default
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop
      
    end select
!EOC

  end subroutine vmix_put_global_params_real

end module vmix_put_get

