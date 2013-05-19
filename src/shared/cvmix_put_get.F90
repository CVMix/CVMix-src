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
                                     cvmix_global_params_type,  &
                                     cvmix_bkgnd_params_type,   &
                                     cvmix_conv_params_type,    &
                                     cvmix_ddiff_params_type,   &
                                     cvmix_shear_params_type
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
    module procedure cvmix_put_bkgnd_real    ! untested
    module procedure cvmix_put_bkgnd_real_1D
    module procedure cvmix_put_bkgnd_real_2D ! untested
    module procedure cvmix_put_conv_real
    module procedure cvmix_put_ddiff_real
    module procedure cvmix_put_shear_real
    module procedure cvmix_put_shear_str
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
        
      case ('diff')
      if (.not.associated(CVmix_vars%diff_iface)) then
        allocate(CVmix_vars%diff_iface(nlev+1,2))
      end if
      if (present(opts)) then
        if (trim(opts).eq.'col2') then
          CVmix_vars%diff_iface(:,2) = val
        else
          CVmix_vars%diff_iface(:,1) = val
        end if
      end if
      
      case ('visc')
      if (.not.associated(CVmix_vars%visc_iface)) then
        allocate(CVmix_vars%visc_iface(nlev+1))
      end if
      CVmix_vars%visc_iface(:) = val

      case default
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
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
    real(cvmix_r8),              intent(in) :: val
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
      case ('diff')
      if (.not.associated(CVmix_vars%diff_iface)) then
        allocate(CVmix_vars%diff_iface(nlev+1,2))
      end if
      if (present(opts)) then
        if (trim(opts).eq.'col2') then
          CVmix_vars%diff_iface(:,2) = val
        else
          CVmix_vars%diff_iface(:,1) = val
        end if
      end if
      
      case ('visc')
      if (.not.associated(CVmix_vars%visc_iface)) then
        allocate(CVmix_vars%visc_iface(nlev+1))
      end if
      CVmix_vars%visc_iface(:) = val

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
        if (trim(opts).eq.'col2') then
          CVmix_vars%diff_iface(:,2) = val
        else
          CVmix_vars%diff_iface(:,1) = val
        end if
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

      case ('zw')
      if (.not.associated(CVmix_vars%z_iface)) then
        allocate(CVmix_vars%z_iface(nlev+1))
      end if
      CVmix_vars%z_iface(:) = val

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

! !IROUTINE: cvmix_put_bkgnd_real
! !INTERFACE:

  subroutine cvmix_put_bkgnd_real(CVmix_bkgnd_params, varname, val)

! !DESCRIPTION:
!  Write a real value into a cvmix\_bkgnd\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    real(cvmix_r8),   intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_bkgnd_params_type), intent(inout) :: CVmix_bkgnd_params
!EOP
!BOC

    select case (trim(varname))
      case ('static_visc')
        if (.not.allocated(CVmix_bkgnd_params%static_visc)) then
          allocate(CVmix_bkgnd_params%static_visc(1,1))
          CVmix_bkgnd_params%lvary_horizontal=.false.
          CVmix_bkgnd_params%lvary_vertical=.false.
        else
          print*, "WARNING: overwriting static_visc in cvmix_bkgnd_params_type."
        end if
        CVmix_bkgnd_params%static_visc(:,:) = val

      case ('static_diff')
        if (.not.allocated(CVmix_bkgnd_params%static_diff)) then
          allocate(CVmix_bkgnd_params%static_diff(1,1))
          CVmix_bkgnd_params%lvary_horizontal=.false.
          CVmix_bkgnd_params%lvary_vertical=.false.
        else
          print*, "WARNING: overwriting static_diff in cvmix_bkgnd_params_type."
        end if
        CVmix_bkgnd_params%static_diff(:,:) = val

      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end subroutine cvmix_put_bkgnd_real

!BOP

! !IROUTINE: cvmix_put_bkgnd_real_1D
! !INTERFACE:

  subroutine cvmix_put_bkgnd_real_1D(CVmix_bkgnd_params, varname, val, &
                                    ncol, nlev)

! !DESCRIPTION:
!  Write an array of real values into a cvmix\_bkgnd\_params\_type variable.
!  You must use \verb|opt='horiz'| to specify that the field varies in the
!  horizontal direction, otherwise it is assumed to vary in the vertical.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*),             intent(in) :: varname
    real(cvmix_r8), dimension(:), intent(in) :: val
    integer, optional,            intent(in) :: ncol, nlev

! !OUTPUT PARAMETERS:
    type (cvmix_bkgnd_params_type), intent(inout) :: CVmix_bkgnd_params
!EOP
!BOC

    ! Local vars
    integer, dimension(2) :: dims
    integer               :: data_dims
    logical               :: lvary_horizontal

    ! Error checking to make sure dimension is specified
    if ((.not.present(ncol)).and.(.not.present(nlev))) then
      print*, "ERROR: when putting 1D data in cvmix_bkgnd_params_type ", &
              "you must specify nlev or ncol!"
      stop 1
    end if

    if ((present(ncol)).and.(present(nlev))) then
      print*, "ERROR: when putting 1D data in cvmix_bkgnd_params_type ", &
              "you can not specify both nlev or ncol!"
      stop 1
    end if

    data_dims = size(val)
    if (present(ncol)) then
      if (data_dims.gt.ncol) then
        print*, "ERROR: data array is bigger than number of columns specified."
        stop 1
      end if
      lvary_horizontal=.true.
      dims(1) = ncol
      dims(2) = 1
    else
      if (data_dims.gt.nlev+1) then
        print*, "ERROR: data array is bigger than number of levels specified."
        stop 1
      end if
      lvary_horizontal=.false.
      dims(1) = 1
      dims(2) = nlev+1
    end if

    select case (trim(varname))
      case ('static_visc')
        if (.not.allocated(CVmix_bkgnd_params%static_visc)) then
          allocate(CVmix_bkgnd_params%static_visc(dims(1),dims(2)))
          CVmix_bkgnd_params%lvary_horizontal = lvary_horizontal
          CVmix_bkgnd_params%lvary_vertical = .not.lvary_horizontal
        else
          print*, "WARNING: overwriting static_visc in cvmix_bkgnd_params_type."
        end if
        if (any(shape(CVmix_bkgnd_params%static_visc).ne.dims)) then
          print*, "ERROR: dimensions of static_visc do not match what was ", &
                  "sent to cvmix_put"
          stop 1
        end if
        if (lvary_horizontal) then
          CVmix_bkgnd_params%static_visc(:,1)           = 0_cvmix_r8
          CVmix_bkgnd_params%static_visc(1:data_dims,1) = val
        else
          CVmix_bkgnd_params%static_visc(1,:)           = 0_cvmix_r8
          CVmix_bkgnd_params%static_visc(1,1:data_dims) = val
        end if

      case ('static_diff')
        if (.not.allocated(CVmix_bkgnd_params%static_diff)) then
          allocate(CVmix_bkgnd_params%static_diff(dims(1),dims(2)))
          CVmix_bkgnd_params%lvary_horizontal = lvary_horizontal
          CVmix_bkgnd_params%lvary_vertical = .not.lvary_horizontal
        else
          print*, "WARNING: overwriting static_diff in cvmix_bkgnd_params_type."
        end if
        if (any(shape(CVmix_bkgnd_params%static_diff).ne.dims)) then
          print*, "ERROR: dimensions of static_diff do not match what was ", &
                  "sent to cvmix_put"
          stop 1
        end if
        if (lvary_horizontal) then
          CVmix_bkgnd_params%static_diff(:,1)           = 0_cvmix_r8
          CVmix_bkgnd_params%static_diff(1:data_dims,1) = val
        else
          CVmix_bkgnd_params%static_diff(1,:)           = 0_cvmix_r8
          CVmix_bkgnd_params%static_diff(1,1:data_dims) = val
        end if

      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end subroutine cvmix_put_bkgnd_real_1D

!BOP

! !IROUTINE: cvmix_put_bkgnd_real_2D
! !INTERFACE:

  subroutine cvmix_put_bkgnd_real_2D(CVmix_bkgnd_params, varname, val, &
                                    ncol, nlev)

! !DESCRIPTION:
!  Write a 2D array of real values into a cvmix\_bkgnd\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*),               intent(in) :: varname
    real(cvmix_r8), dimension(:,:), intent(in) :: val
    integer,                        intent(in) :: ncol, nlev

! !OUTPUT PARAMETERS:
    type (cvmix_bkgnd_params_type), intent(out) :: CVmix_bkgnd_params
!EOP
!BOC

    ! Local vars
    integer, dimension(2) :: dims, data_dims

    dims      = (/ncol, nlev+1/)
    data_dims = shape(val)

    if (any(data_dims.gt.dims)) then
      print*, "ERROR: data being put in cvmix_bkgnd_params_type is larger ", &
              "than (ncol, nlev+1)"
      stop 1
    end if

    select case (trim(varname))
      case ('static_visc')
        if (.not.allocated(CVmix_bkgnd_params%static_visc)) then
          allocate(CVmix_bkgnd_params%static_visc(dims(1),dims(2)))
          CVmix_bkgnd_params%lvary_horizontal=.true.
          CVmix_bkgnd_params%lvary_vertical=.true.
        else
          print*, "WARNING: overwriting static_visc in cvmix_bkgnd_params_type."
        end if
        if (any(shape(CVmix_bkgnd_params%static_visc).ne.dims)) then
          print*, "ERROR: dimensions of static_visc do not match what was ", &
                  "sent to cvmix_put"
          stop 1
        end if
        CVmix_bkgnd_params%static_visc = 0.0_cvmix_r8
        CVmix_bkgnd_params%static_visc(1:data_dims(1), 1:data_dims(2)) = val

      case ('static_diff')
        if (.not.allocated(CVmix_bkgnd_params%static_diff)) then
          allocate(CVmix_bkgnd_params%static_diff(dims(1),dims(2)))
          CVmix_bkgnd_params%lvary_horizontal=.true.
          CVmix_bkgnd_params%lvary_vertical=.true.
        else
          print*, "WARNING: overwriting static_diff in cvmix_bkgnd_params_type."
        end if
        if (any(shape(CVmix_bkgnd_params%static_diff).ne.dims)) then
          print*, "ERROR: dimensions of static_diff do not match what was ", &
                  "sent to cvmix_put"
          stop 1
        end if
        CVmix_bkgnd_params%static_diff = 0.0_cvmix_r8
        CVmix_bkgnd_params%static_diff(1:data_dims(1), 1:data_dims(2)) = val

      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end subroutine cvmix_put_bkgnd_real_2D

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

! !IROUTINE: cvmix_put_ddiff_real
! !INTERFACE:

  subroutine cvmix_put_ddiff_real(CVmix_ddiff_params, varname, val)

! !DESCRIPTION:
!  Write a real value into a cvmix\_ddiff\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    real(cvmix_r8),   intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_ddiff_params_type), intent(inout) :: CVmix_ddiff_params
!EOP
!BOC

    select case (trim(varname))
      case ('strat_param_max')
        CVmix_ddiff_params%strat_param_max = val
      case ('ddiff_exp1')
        CVmix_ddiff_params%ddiff_exp1 = val
      case ('ddiff_exp2')
        CVmix_ddiff_params%ddiff_exp2 = val
      case ('kappa_ddiff_param1')
        CVmix_ddiff_params%kappa_ddiff_param1 = val
      case ('kappa_ddiff_param2')
        CVmix_ddiff_params%kappa_ddiff_param2 = val
      case ('kappa_ddiff_param3')
        CVmix_ddiff_params%kappa_ddiff_param3 = val
      case ('kappa_ddiff_t')
        CVmix_ddiff_params%kappa_ddiff_t = val
      case ('kappa_ddiff_s')
        CVmix_ddiff_params%kappa_ddiff_s = val
      case ('mol_diff')
        CVmix_ddiff_params%mol_diff = val
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end subroutine cvmix_put_ddiff_real

!BOP

! !IROUTINE: cvmix_put_shear_real
! !INTERFACE:

  subroutine cvmix_put_shear_real(CVmix_shear_params, varname, val)

! !DESCRIPTION:
!  Write a real value into a cvmix\_shear\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    real(cvmix_r8),   intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_shear_params_type), intent(inout) :: CVmix_shear_params
!EOP
!BOC

    select case (trim(varname))
      case ('PP_nu_zero')
        CVmix_shear_params%PP_nu_zero = val
      case ('PP_alpha')
        CVmix_shear_params%PP_alpha = val
      case ('PP_exp')
        CVmix_shear_params%PP_exp = val
      case ('KPP_nu_zero')
        CVmix_shear_params%KPP_nu_zero = val
      case ('KPP_Ri_zero')
        CVmix_shear_params%KPP_Ri_zero = val
      case ('KPP_exp')
        CVmix_shear_params%KPP_exp = val
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end subroutine cvmix_put_shear_real

!BOP

! !IROUTINE: cvmix_put_shear_str
! !INTERFACE:

  subroutine cvmix_put_shear_str(CVmix_shear_params, varname, val)

! !DESCRIPTION:
!  Write a string into a cvmix\_shear\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_shear_params_type), intent(inout) :: CVmix_shear_params
!EOP
!BOC

    select case (trim(varname))
      case ('mix_scheme')
        CVmix_shear_params%mix_scheme = val
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end subroutine cvmix_put_shear_str

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
        
      case default
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select
!EOC

  end subroutine cvmix_put_global_params_real

end module cvmix_put_get

