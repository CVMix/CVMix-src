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

  subroutine cvmix_put_int(CVmix_vars, varname, val)

! !DESCRIPTION:
!  Write an integer value into a cvmix\_data\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    integer,          intent(in) :: val

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
      case default
        ! All other scalars are real(cvmix_r8)
        call cvmix_put_real(CVmix_vars, varname, real(val,cvmix_r8))
    end select
!EOC

  end subroutine cvmix_put_int

!BOP

! !IROUTINE: cvmix_put_real
! !INTERFACE:

  subroutine cvmix_put_real(CVmix_vars, varname, val)

! !DESCRIPTION:
!  Write a real value into a cvmix\_data\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*),           intent(in) :: varname
    real(cvmix_r8),             intent(in) :: val

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
      case ('OceanDepth','ocn_depth','depth')
        CVmix_vars%OceanDepth = val
      case ('BoundaryLayerDepth','OBL_depth')
        CVmix_vars%BoundaryLayerDepth = val
      case ('kOBL_depth')
        CVmix_vars%kOBL_depth = val
      case ('SurfaceHeight','SeaSurfaceHeight','surf_hgt')
        CVmix_vars%SeaSurfaceHeight = val
      case ('SurfaceFriction','surf_fric')
        CVmix_vars%SurfaceFriction = val
      case ('SurfaceBuoyancy','SurfaceBuoyancyForcing','surf_buoy')
        CVmix_vars%SurfaceBuoyancyForcing = val
      case ('lat','latitude')
        CVmix_vars%lat = val
      case ('lon','longitude')
        CVmix_vars%lon = val
      case ('Coriolis')
        CVmix_vars%Coriolis = val

      case ('Mdiff')
      if (.not.associated(CVmix_vars%Mdiff_iface)) then
        allocate(CVmix_vars%Mdiff_iface(nlev+1))
      end if
      CVmix_vars%Mdiff_iface(:) = val

      case ('Tdiff')
      if (.not.associated(CVmix_vars%Tdiff_iface)) then
        allocate(CVmix_vars%Tdiff_iface(nlev+1))
      end if
      CVmix_vars%Tdiff_iface(:) = val

      case ('Sdiff')
      if (.not.associated(CVmix_vars%Sdiff_iface)) then
        allocate(CVmix_vars%Sdiff_iface(nlev+1))
      end if
      CVmix_vars%Sdiff_iface(:) = val

      case ('WaterDensity','dens')
      print*, "WARNING: you are setting the density in all levels to a",      &
              "constant value"
      if (.not.associated(CVmix_vars%WaterDensity_cntr)) then
        allocate(CVmix_vars%WaterDensity_cntr(nlev))
      end if
      CVmix_vars%WaterDensity_cntr(:) = val

      case ('dens_lwr')
      print*, "WARNING: you are setting the adiabatic density in all levels", &
              "to a constant value"
      if (.not.associated(CVmix_vars%AdiabWaterDensity_cntr)) then
        allocate(CVmix_vars%AdiabWaterDensity_cntr(nlev))
      end if
      CVmix_vars%AdiabWaterDensity_cntr(:) = val

      case ('Richardson','ShearRichardson','RichardsonNumber',                &
            'ShearRichardsonNumer','Ri','Ri_iface')
      print*, "WARNING: you are setting the Richardson number in all", &
              "levels to a constant value"
      if (.not.associated(CVmix_vars%ShearRichardson_iface)) then
        allocate(CVmix_vars%ShearRichardson_iface(nlev+1))
      end if
      CVmix_vars%ShearRichardson_iface(:) = val

      case ('BulkRichardson','BulkRichardsonNumber','Rib','Ri_bulk')
      print*, "WARNING: you are setting the bulk Richardson number in all",  &
              "levels to a constant value"
      if (.not.associated(CVmix_vars%BulkRichardson_cntr)) then
        allocate(CVmix_vars%BulkRichardson_cntr(nlev))
      end if
      CVmix_vars%BulkRichardson_cntr(:) = val

      case ('dz','dzt')
      print*, "WARNING: you are setting the cell thickness in all levels to", &
              "a constant value"
      if (.not.associated(CVmix_vars%dzt)) then
        allocate(CVmix_vars%dzt(nlev))
      end if
      CVmix_vars%dzt(:) = val

      case ('dzw', 'dzw_iface')
      print*, "WARNING: you are setting the cell midpoint to midpoint", &
              "distance in all levels to a constant value"
      if (.not.associated(CVmix_vars%dzw)) then
        allocate(CVmix_vars%dzw(nlev+1))
      end if
      CVmix_vars%dzw(:) = val

      case ('Buoyancy','SqrBuoyancy','SqrBuoyancyFreq','buoy', 'buoy_iface')
      print*, "WARNING: you are setting the buoyancy in all levels to a", &
              "constant value"
      if (.not.associated(CVmix_vars%SqrBuoyancyFreq_iface)) then
        allocate(CVmix_vars%SqrBuoyancyFreq_iface(nlev+1))
      end if
      CVmix_vars%SqrBuoyancyFreq_iface(:) = val

      case ('kpp_transport','kpp_nonlocal','nonlocal_transport')
      if (.not.associated(CVmix_vars%kpp_Tnonlocal_iface)) then
        allocate(CVmix_vars%kpp_Tnonlocal_iface(nlev+1))
      end if
      if (.not.associated(CVmix_vars%kpp_Snonlocal_iface)) then
        allocate(CVmix_vars%kpp_Snonlocal_iface(nlev+1))
      end if
      CVmix_vars%kpp_Tnonlocal_iface(:) = val
      CVmix_vars%kpp_Snonlocal_iface(:) = val

      case ('kpp_Ttransport','kpp_Tnonlocal','Tnonlocal')
      if (.not.associated(CVmix_vars%kpp_Tnonlocal_iface)) then
        allocate(CVmix_vars%kpp_Tnonlocal_iface(nlev+1))
      end if
      CVmix_vars%kpp_Tnonlocal_iface(:) = val

      case ('kpp_Stransport','kpp_Snonlocal','Snonlocal')
      if (.not.associated(CVmix_vars%kpp_Snonlocal_iface)) then
        allocate(CVmix_vars%kpp_Snonlocal_iface(nlev+1))
      end if
      CVmix_vars%kpp_Snonlocal_iface(:) = val

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

  subroutine cvmix_put_real_1D(CVmix_vars, varname, val)

! !DESCRIPTION:
!  Write an array of real values into a cvmix\_data\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*),             intent(in) :: varname
    real(cvmix_r8), dimension(:), intent(in) :: val

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
      case ('Mdiff')
      if (.not.associated(CVmix_vars%Mdiff_iface)) then
        allocate(CVmix_vars%Mdiff_iface(nlev+1))
      end if
      CVmix_vars%Mdiff_iface(:) = val

      case ('Tdiff')
      if (.not.associated(CVmix_vars%Tdiff_iface)) then
        allocate(CVmix_vars%Tdiff_iface(nlev+1))
      end if
      CVmix_vars%Tdiff_iface(:) = val

      case ('Sdiff')
      if (.not.associated(CVmix_vars%Sdiff_iface)) then
        allocate(CVmix_vars%Sdiff_iface(nlev+1))
      end if
      CVmix_vars%Sdiff_iface(:) = val

      case ('WaterDensity','dens')
      if (.not.associated(CVmix_vars%WaterDensity_cntr)) then
        allocate(CVmix_vars%WaterDensity_cntr(nlev))
      end if
      CVmix_vars%WaterDensity_cntr(:) = val

      case ('dens_lwr')
      if (.not.associated(CVmix_vars%AdiabWaterDensity_cntr)) then
        allocate(CVmix_vars%AdiabWaterDensity_cntr(nlev))
      end if
      CVmix_vars%AdiabWaterDensity_cntr(:) = val

      case ('Richardson','ShearRichardson','RichardsonNumber',                &
            'ShearRichardsonNumer','Ri','Ri_iface')
      if (.not.associated(CVmix_vars%ShearRichardson_iface)) then
        allocate(CVmix_vars%ShearRichardson_iface(nlev+1))
      end if
      CVmix_vars%ShearRichardson_iface(:) = val

      case ('BulkRichardson','BulkRichardsonNumber','Rib','Ri_bulk')
      if (.not.associated(CVmix_vars%BulkRichardson_cntr)) then
        allocate(CVmix_vars%BulkRichardson_cntr(nlev))
      end if
      CVmix_vars%BulkRichardson_cntr(:) = val

      case ('z','zt','zt_cntr')
      if (.not.associated(CVmix_vars%zt_cntr)) then
        allocate(CVmix_vars%zt_cntr(nlev))
      end if
      CVmix_vars%zt_cntr(:) = val

      case ('dz','dzt')
      if (.not.associated(CVmix_vars%dzt)) then
        allocate(CVmix_vars%dzt(nlev))
      end if
      CVmix_vars%dzt(:) = val

      case ('zw','zw_iface')
      if (.not.associated(CVmix_vars%zw_iface)) then
        allocate(CVmix_vars%zw_iface(nlev+1))
      end if
      CVmix_vars%zw_iface(:) = val

      case ('dzw')
      if (.not.associated(CVmix_vars%dzw)) then
        allocate(CVmix_vars%dzw(nlev+1))
      end if
      CVmix_vars%dzw(:) = val

      case ('Buoyancy','SqrBuoyancy','SqrBuoyancyFreq','buoy', 'buoy_iface')
      if (.not.associated(CVmix_vars%SqrBuoyancyFreq_iface)) then
        allocate(CVmix_vars%SqrBuoyancyFreq_iface(nlev+1))
      end if
      CVmix_vars%SqrBuoyancyFreq_iface(:) = val

      case ('kpp_transport','kpp_nonlocal','nonlocal_transport')
      if (.not.associated(CVmix_vars%kpp_Tnonlocal_iface)) then
        allocate(CVmix_vars%kpp_Tnonlocal_iface(nlev+1))
      end if
      if (.not.associated(CVmix_vars%kpp_Snonlocal_iface)) then
        allocate(CVmix_vars%kpp_Snonlocal_iface(nlev+1))
      end if
      CVmix_vars%kpp_Tnonlocal_iface(:) = val
      CVmix_vars%kpp_Snonlocal_iface(:) = val

      case ('kpp_Ttransport','kpp_Tnonlocal','Tnonlocal')
      if (.not.associated(CVmix_vars%kpp_Tnonlocal_iface)) then
        allocate(CVmix_vars%kpp_Tnonlocal_iface(nlev+1))
      end if
      CVmix_vars%kpp_Tnonlocal_iface(:) = val

      case ('kpp_Stransport','kpp_Snonlocal','Snonlocal')
      if (.not.associated(CVmix_vars%kpp_Snonlocal_iface)) then
        allocate(CVmix_vars%kpp_Snonlocal_iface(nlev+1))
      end if
      CVmix_vars%kpp_Snonlocal_iface(:) = val

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
      case ('fw_rho','FreshWaterDensity')
        CVmix_params%FreshWaterDensity = val
      case ('sw_rho','SaltWaterDensity')
        CVmix_params%SaltWaterDensity= val
        
      case default
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select
!EOC

  end subroutine cvmix_put_global_params_real

end module cvmix_put_get

