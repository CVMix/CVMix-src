module cvmix_put_get

!BOP
!\newpage
! !MODULE: cvmix_put_get
!
! !DESCRIPTION:
!  This module contains routines to pack data into the cvmix datatypes
!  (allocating memory as necessary) and then unpack the data out. If we switch
!  to pointers, the pack will just point at the right target and the unpack
!  will be un-necessary.
!\\
!\\

! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !USES:

   use cvmix_kinds_and_types, only : cvmix_r8,                  &
                                     cvmix_data_type,           &
                                     cvmix_global_params_type
!EOP

  implicit none
  private
  save

!BOP

! !PUBLIC MEMBER FUNCTIONS:
  public :: cvmix_put

  interface cvmix_put
    module procedure cvmix_put_int
    module procedure cvmix_put_real
    module procedure cvmix_put_real_1D
    module procedure cvmix_put_global_params_int
    module procedure cvmix_put_global_params_real
  end interface cvmix_put
!EOP

contains

!BOP

! !IROUTINE: cvmix_put_int
! !INTERFACE:

  subroutine cvmix_put_int(CVmix_vars, varname, val, opts)

! !DESCRIPTION:
!  Write an integer value into a cvmix\_data\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*),           intent(in) :: varname
    integer,                    intent(in) :: val
    character(len=*), optional, intent(in) :: opts

! !OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars
!EOP
!BOC

    ! Local variables
    integer :: nlev
  
    nlev = CVmix_vars%nlev

    if ((trim(varname).ne.'nlev').and.(nlev.eq.-1)) then
      print*, "ERROR: you must specify the number of levels before ", &
              "you can pack data into a cvmix_data_type!"
      print*, "You tried to set ", trim(varname)
      stop 1
    end if
    
    select case (trim(varname))
      case ('nlev')
        CVmix_vars%nlev = val
      case ('kOBL_depth')
        CVmix_vars%kOBL_depth = val
      case default
        ! All other scalars are real(cvmix_r8)
        call cvmix_put_real(CVmix_vars, varname, real(val,cvmix_r8), opts)
    end select
!EOC

  end subroutine cvmix_put_int

!BOP

! !IROUTINE: cvmix_put_real
! !INTERFACE:

  subroutine cvmix_put_real(CVmix_vars, varname, val, opts)

! !DESCRIPTION:
!  Write a real value into a cvmix\_data\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*),           intent(in) :: varname
    real(cvmix_r8),             intent(in) :: val
    character(len=*), optional, intent(in) :: opts

! !OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars
!EOP
!BOC

    ! Local variables
    integer :: nlev
  
    nlev = CVmix_vars%nlev

    if (nlev.eq.-1) then
      print*, "ERROR: you must specify the number of levels before ", &
              "you can pack data into a cvmix_data_type!"
      print*, "You tried to set ", trim(varname)
      stop 1
    end if
    
    select case (trim(varname))
      case ('ocn_depth','depth')
        CVmix_vars%ocn_depth = val
      case ('OBL_depth')
        CVmix_vars%OBL_depth = val
      case ('surf_hgt')
        CVmix_vars%surf_hgt = val
      case ('surf_fric')
        CVmix_vars%surf_fric = val
      case ('surf_buoy')
        CVmix_vars%surf_buoy = val
      case ('lat')
        CVmix_vars%lat = val
      case ('lon')
        CVmix_vars%lon = val
      case ('Coriolis')
        CVmix_vars%Coriolis = val
      case ('diff')
        if (.not.associated(CVmix_vars%diff_iface)) then
          allocate(CVmix_vars%diff_iface(nlev+1,2))
        end if
        if (present(opts)) then
          select case (trim(opts))
            case ('col1')
              CVmix_vars%diff_iface(:,1) = val
            case ('col2')
              CVmix_vars%diff_iface(:,2) = val
            case DEFAULT
              print*, "WARNING: ignoring opts = ", trim(opts)
              CVmix_vars%diff_iface      = val
          end select
        end if
      case ('visc')
      if (.not.associated(CVmix_vars%visc_iface)) then
        allocate(CVmix_vars%visc_iface(nlev+1))
      end if
      CVmix_vars%visc_iface(:) = val

      case ('dens')
      print*, "WARNING: you are setting the density in all levels to a", &
              "constant value"
      if (.not.associated(CVmix_vars%dens)) then
        allocate(CVmix_vars%dens(nlev))
      end if
      CVmix_vars%dens(:) = val

      case ('dens_lwr')
      print*, "WARNING: you are setting the density in all levels to a", &
              "constant value"
      if (.not.associated(CVmix_vars%dens_lwr)) then
        allocate(CVmix_vars%dens_lwr(nlev))
      end if
      CVmix_vars%dens_lwr(:) = val

      case ('Ri', 'Ri_iface')
      print*, "WARNING: you are setting the Richardson number in all", &
              "levels to a constant value"
      if (.not.associated(CVmix_vars%Ri_iface)) then
        allocate(CVmix_vars%Ri_iface(nlev+1))
      end if
      CVmix_vars%Ri_iface(:) = val

      case ('dz', 'dzt')
      print*, "WARNING: you are setting the cell thickness in all levels to", &
              "a constant value"
      if (.not.associated(CVmix_vars%dzt)) then
        allocate(CVmix_vars%dzt(nlev))
      end if
      CVmix_vars%dzt(:) = val

      case ('dzw', 'dzw_iface')
      print*, "WARNING: you are setting the cell midpoint to midpoint", &
              "distance in all levels to a constant value"
      if (.not.associated(CVmix_vars%dzw_iface)) then
        allocate(CVmix_vars%dzw_iface(nlev+1))
      end if
      CVmix_vars%dzw_iface(:) = val

      case ('buoy', 'buoy_iface')
      print*, "WARNING: you are setting the buoyancy in all levels to a", &
              "constant value"
      if (.not.associated(CVmix_vars%buoy_iface)) then
        allocate(CVmix_vars%buoy_iface(nlev+1))
      end if
      CVmix_vars%buoy_iface(:) = val

      case ('strat_param_num')
      print*, "WARNING: you are setting the numerator of the statification", & 
              "parameter in all levels to a constant value"
      if (.not.associated(CVmix_vars%strat_param_num)) then
        allocate(CVmix_vars%strat_param_num(nlev))
      end if
      CVmix_vars%strat_param_num(:) = val

      case ('strat_param_denom')
      print*, "WARNING: you are setting the denominator of the statification",& 
              "parameter in all levels to a constant value"
      if (.not.associated(CVmix_vars%strat_param_denom)) then
        allocate(CVmix_vars%strat_param_denom(nlev))
      end if
      CVmix_vars%strat_param_denom(:) = val

      case default
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select
!EOC

  end subroutine cvmix_put_real

!BOP

! !IROUTINE: cvmix_put_real_1D
! !INTERFACE:

  subroutine cvmix_put_real_1D(CVmix_vars, varname, val, opts)

! !DESCRIPTION:
!  Write an array of real values into a cvmix\_data\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*),             intent(in) :: varname
    real(cvmix_r8), dimension(:), intent(in) :: val
    character(len=*), optional,   intent(in) :: opts

! !OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars
!EOP
!BOC

    ! Local variables
    integer :: nlev
  
    nlev = CVmix_vars%nlev

    if (nlev.eq.-1) then
      print*, "ERROR: you must specify the number of levels before ", &
              "you can pack data into a cvmix_data_type!"
      print*, "You tried to set ", trim(varname)
      stop 1
    end if
    
    select case (trim(varname))
      case ('diff')
      if (.not.associated(CVmix_vars%diff_iface)) then
        allocate(CVmix_vars%diff_iface(nlev+1,2))
      end if
      if (present(opts)) then
        select case (trim(opts))
          case ('col1')
            CVmix_vars%diff_iface(:,1) = val
          case ('col2')
            CVmix_vars%diff_iface(:,2) = val
          case DEFAULT
            print*, "WARNING: ignoring opts = ", trim(opts)
            CVmix_vars%diff_iface(:,1) = val
            CVmix_vars%diff_iface(:,2) = val
        end select
      end if

      case ('visc')
      if (.not.associated(CVmix_vars%visc_iface)) then
        allocate(CVmix_vars%visc_iface(nlev+1))
      end if
      CVmix_vars%visc_iface(:) = val

      case ('dens')
      if (.not.associated(CVmix_vars%dens)) then
        allocate(CVmix_vars%dens(nlev))
      end if
      CVmix_vars%dens(:) = val

      case ('dens_lwr')
      if (.not.associated(CVmix_vars%dens_lwr)) then
        allocate(CVmix_vars%dens_lwr(nlev))
      end if
      CVmix_vars%dens_lwr(:) = val

      case ('Ri', 'Ri_iface')
      if (.not.associated(CVmix_vars%Ri_iface)) then
        allocate(CVmix_vars%Ri_iface(nlev+1))
      end if
      CVmix_vars%Ri_iface(:) = val

      case ('z', 'zt')
      if (.not.associated(CVmix_vars%zt)) then
        allocate(CVmix_vars%zt(nlev))
      end if
      CVmix_vars%zt(:) = val

      case ('dz', 'dzt')
      if (.not.associated(CVmix_vars%dzt)) then
        allocate(CVmix_vars%dzt(nlev))
      end if
      CVmix_vars%dzt(:) = val

      case ('zw', 'zw_iface')
      if (.not.associated(CVmix_vars%zw_iface)) then
        allocate(CVmix_vars%zw_iface(nlev+1))
      end if
      CVmix_vars%zw_iface(:) = val

      case ('dzw', 'dzw_iface')
      if (.not.associated(CVmix_vars%dzw_iface)) then
        allocate(CVmix_vars%dzw_iface(nlev+1))
      end if
      CVmix_vars%dzw_iface(:) = val

      case ('buoy', 'buoy_iface')
      if (.not.associated(CVmix_vars%buoy_iface)) then
        allocate(CVmix_vars%buoy_iface(nlev+1))
      end if
      CVmix_vars%buoy_iface(:) = val

      case ('strat_param_num')
      if (.not.associated(CVmix_vars%strat_param_num)) then
        allocate(CVmix_vars%strat_param_num(nlev))
      end if
      CVmix_vars%strat_param_num(:) = val

      case ('strat_param_denom')
      if (.not.associated(CVmix_vars%strat_param_denom)) then
        allocate(CVmix_vars%strat_param_denom(nlev))
      end if
      CVmix_vars%strat_param_denom(:) = val

      case default
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end subroutine cvmix_put_real_1D

!BOP

! !IROUTINE: cvmix_put_global_params_int
! !INTERFACE:

  subroutine cvmix_put_global_params_int(CVmix_params, varname, val)

! !DESCRIPTION:
!  Write an integer value into a cvmix\_global\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    integer,          intent(in) :: val

! !OUTPUT PARAMETERS:
    type (cvmix_global_params_type), intent(inout) :: CVmix_params
!EOP
!BOC

    select case (trim(varname))
      case ('max_nlev')
        CVmix_params%max_nlev = val
        
      case default
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select
!EOC

  end subroutine cvmix_put_global_params_int

!BOP

! !IROUTINE: cvmix_put_global_params_real
! !INTERFACE:

  subroutine cvmix_put_global_params_real(CVmix_params, varname, val)

! !DESCRIPTION:
!  Write a real value into a cvmix\_global\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    real(cvmix_r8),   intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_global_params_type), intent(inout) :: CVmix_params
!EOP
!BOC

    select case (trim(varname))
      case ('prandtl')
        CVmix_params%prandtl = val
      case ('fw_rho')
        CVmix_params%fw_rho  = val
      case ('sw_rho')
        CVmix_params%sw_rho  = val
        
      case default
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select
!EOC

  end subroutine cvmix_put_global_params_real

end module cvmix_put_get

