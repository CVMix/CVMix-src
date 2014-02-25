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
                                     cvmix_strlen,              &
                                     cvmix_data_type,           &
                                     cvmix_global_params_type
!EOP

  implicit none
  private
  save

!BOP

! !PUBLIC MEMBER FUNCTIONS:
  public :: cvmix_put
  public :: cvmix_att_name

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
    
    select case (trim(cvmix_att_name(varname)))
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
    
    select case (trim(cvmix_att_name(varname)))
      case ('OceanDepth')
        CVmix_vars%OceanDepth = val
      case ('BoundaryLayerDepth')
        CVmix_vars%BoundaryLayerDepth = val
      case ('SeaSurfaceHeight')
        CVmix_vars%SeaSurfaceHeight = val
      case ('SurfaceFriction')
        CVmix_vars%SurfaceFriction = val
      case ("SurfaceBuoyancyForcing")
        CVmix_vars%SurfaceBuoyancyForcing = val
      case ("lat")
        CVmix_vars%lat = val
      case ("lon")
        CVmix_vars%lon = val
      case ("Coriolis")
        CVmix_vars%Coriolis = val
      case ("kOBL_depth")
        CVmix_vars%kOBL_depth = val

      case ("dzw")
        print*, "WARNING: you are setting the cell midpoint to midpoint",     &
                "distance in all levels to a constant value"
        if (.not.associated(CVmix_vars%dzw)) then
          allocate(CVmix_vars%dzw(nlev+1))
        end if
        CVmix_vars%dzw(:) = val
      case ("Mdiff_iface")
        if (.not.associated(CVmix_vars%Mdiff_iface)) then
          allocate(CVmix_vars%Mdiff_iface(nlev+1))
        end if
        CVmix_vars%Mdiff_iface(:) = val
      case ("Tdiff_iface")
        if (.not.associated(CVmix_vars%Tdiff_iface)) then
          allocate(CVmix_vars%Tdiff_iface(nlev+1))
        end if
        CVmix_vars%Tdiff_iface(:) = val
      case ("Sdiff_iface")
        if (.not.associated(CVmix_vars%Sdiff_iface)) then
          allocate(CVmix_vars%Sdiff_iface(nlev+1))
        end if
        CVmix_vars%Sdiff_iface(:) = val
      case ("ShearRichardson_iface")
        print*, "WARNING: you are setting the Richardson number in all",      &
                "levels to a constant value"
        if (.not.associated(CVmix_vars%ShearRichardson_iface)) then
          allocate(CVmix_vars%ShearRichardson_iface(nlev+1))
        end if
        CVmix_vars%ShearRichardson_iface(:) = val
      case ("SqrBuoyancyFreq_iface")
        print*, "WARNING: you are setting the buoyancy in all levels to a", &
                "constant value"
        if (.not.associated(CVmix_vars%SqrBuoyancyFreq_iface)) then
          allocate(CVmix_vars%SqrBuoyancyFreq_iface(nlev+1))
        end if
        CVmix_vars%SqrBuoyancyFreq_iface(:) = val
      case ("kpp_nonlocal_iface")
        if (.not.associated(CVmix_vars%kpp_Tnonlocal_iface)) then
          allocate(CVmix_vars%kpp_Tnonlocal_iface(nlev+1))
        end if
        if (.not.associated(CVmix_vars%kpp_Snonlocal_iface)) then
          allocate(CVmix_vars%kpp_Snonlocal_iface(nlev+1))
        end if
        CVmix_vars%kpp_Tnonlocal_iface(:) = val
        CVmix_vars%kpp_Snonlocal_iface(:) = val
      case ("kpp_Tnonlocal_iface")
        if (.not.associated(CVmix_vars%kpp_Tnonlocal_iface)) then
          allocate(CVmix_vars%kpp_Tnonlocal_iface(nlev+1))
        end if
        CVmix_vars%kpp_Tnonlocal_iface(:) = val
      case ("kpp_Snonlocal_iface")
        if (.not.associated(CVmix_vars%kpp_Snonlocal_iface)) then
          allocate(CVmix_vars%kpp_Snonlocal_iface(nlev+1))
        end if
        CVmix_vars%kpp_Snonlocal_iface(:) = val

      case ("dzt")
        print*, "WARNING: you are setting the cell thickness in all levels",  &
                "to a constant value"
        if (.not.associated(CVmix_vars%dzt)) then
          allocate(CVmix_vars%dzt(nlev))
        end if
        CVmix_vars%dzt(:) = val
      case ("WaterDensity_cntr")
        print*, "WARNING: you are setting the density in all levels to a",    &
                "constant value"
        if (.not.associated(CVmix_vars%WaterDensity_cntr)) then
          allocate(CVmix_vars%WaterDensity_cntr(nlev))
        end if
        CVmix_vars%WaterDensity_cntr(:) = val
      case ("AdiabWaterDensity_cntr")
        print*, "WARNING: you are setting the adiabatic density in all",      &
                "levels to a constant value"
        if (.not.associated(CVmix_vars%AdiabWaterDensity_cntr)) then
          allocate(CVmix_vars%AdiabWaterDensity_cntr(nlev))
        end if
        CVmix_vars%AdiabWaterDensity_cntr(:) = val
      case ("BulkRichardson_cntr")
        print*, "WARNING: you are setting the bulk Richardson number in all", &
                "levels to a constant value"
        if (.not.associated(CVmix_vars%BulkRichardson_cntr)) then
          allocate(CVmix_vars%BulkRichardson_cntr(nlev))
        end if
        CVmix_vars%BulkRichardson_cntr(:) = val
      case ('strat_param_num')
        print*, "WARNING: you are setting the numerator of the",              & 
                "stratification parameter in all levels to a constant value"
        if (.not.associated(CVmix_vars%strat_param_num)) then
          allocate(CVmix_vars%strat_param_num(nlev))
        end if
        CVmix_vars%strat_param_num(:) = val
      case ('strat_param_denom')
        print*, "WARNING: you are setting the denominator of the",            & 
                "stratification parameter in all levels to a constant value"
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
    
    select case (trim(cvmix_att_name(varname)))
      case ("zw_iface")
        if (.not.associated(CVmix_vars%zw_iface)) then
          allocate(CVmix_vars%zw_iface(nlev+1))
        end if
        CVmix_vars%zw_iface(:) = val
      case ("dzw")
        if (.not.associated(CVmix_vars%dzw)) then
          allocate(CVmix_vars%dzw(nlev+1))
        end if
        CVmix_vars%dzw(:) = val
      case ("Mdiff_iface")
        if (.not.associated(CVmix_vars%Mdiff_iface)) then
          allocate(CVmix_vars%Mdiff_iface(nlev+1))
        end if
        CVmix_vars%Mdiff_iface(:) = val
      case ("Tdiff_iface")
        if (.not.associated(CVmix_vars%Tdiff_iface)) then
          allocate(CVmix_vars%Tdiff_iface(nlev+1))
        end if
        CVmix_vars%Tdiff_iface(:) = val
      case ("Sdiff_iface")
        if (.not.associated(CVmix_vars%Sdiff_iface)) then
          allocate(CVmix_vars%Sdiff_iface(nlev+1))
        end if
        CVmix_vars%Sdiff_iface(:) = val
      case ("ShearRichardson_iface")
        if (.not.associated(CVmix_vars%ShearRichardson_iface)) then
          allocate(CVmix_vars%ShearRichardson_iface(nlev+1))
        end if
        CVmix_vars%ShearRichardson_iface(:) = val
      case ("SqrBuoyancyFreq_iface")
        if (.not.associated(CVmix_vars%SqrBuoyancyFreq_iface)) then
          allocate(CVmix_vars%SqrBuoyancyFreq_iface(nlev+1))
        end if
        CVmix_vars%SqrBuoyancyFreq_iface(:) = val
      case ("kpp_nonlocal_iface")
        if (.not.associated(CVmix_vars%kpp_Tnonlocal_iface)) then
          allocate(CVmix_vars%kpp_Tnonlocal_iface(nlev+1))
        end if
        if (.not.associated(CVmix_vars%kpp_Snonlocal_iface)) then
          allocate(CVmix_vars%kpp_Snonlocal_iface(nlev+1))
        end if
        CVmix_vars%kpp_Tnonlocal_iface(:) = val
        CVmix_vars%kpp_Snonlocal_iface(:) = val
      case ("kpp_Tnonlocal_iface")
        if (.not.associated(CVmix_vars%kpp_Tnonlocal_iface)) then
          allocate(CVmix_vars%kpp_Tnonlocal_iface(nlev+1))
        end if
        CVmix_vars%kpp_Tnonlocal_iface(:) = val
      case ("kpp_Snonlocal_iface")
        if (.not.associated(CVmix_vars%kpp_Snonlocal_iface)) then
          allocate(CVmix_vars%kpp_Snonlocal_iface(nlev+1))
        end if
        CVmix_vars%kpp_Snonlocal_iface(:) = val

      case ("zt_cntr")
        if (.not.associated(CVmix_vars%zt_cntr)) then
          allocate(CVmix_vars%zt_cntr(nlev))
        end if
        CVmix_vars%zt_cntr(:) = val
      case ("dzt")
        if (.not.associated(CVmix_vars%dzt)) then
          allocate(CVmix_vars%dzt(nlev))
        end if
        CVmix_vars%dzt(:) = val
      case ("WaterDensity_cntr")
        if (.not.associated(CVmix_vars%WaterDensity_cntr)) then
          allocate(CVmix_vars%WaterDensity_cntr(nlev))
        end if
        CVmix_vars%WaterDensity_cntr(:) = val
      case ("AdiabWaterDensity_cntr")
        if (.not.associated(CVmix_vars%AdiabWaterDensity_cntr)) then
          allocate(CVmix_vars%AdiabWaterDensity_cntr(nlev))
        end if
        CVmix_vars%AdiabWaterDensity_cntr(:) = val
      case ("BulkRichardson_cntr")
        if (.not.associated(CVmix_vars%BulkRichardson_cntr)) then
          allocate(CVmix_vars%BulkRichardson_cntr(nlev))
        end if
        CVmix_vars%BulkRichardson_cntr(:) = val
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
      case ("buoyancy_cntr")
        if (.not.associated(CVmix_vars%buoyancy_cntr)) then
          allocate(CVmix_vars%buoyancy_cntr(nlev))
        end if
        CVmix_vars%buoyancy_cntr(:) = val
      case ("Vx_cntr")
        if (.not.associated(CVmix_vars%Vx_cntr)) then
          allocate(CVmix_vars%Vx_cntr(nlev))
        end if
        CVmix_vars%Vx_cntr(:) = val
      case ("Vy_cntr")
        if (.not.associated(CVmix_vars%Vy_cntr)) then
          allocate(CVmix_vars%Vy_cntr(nlev))
        end if
        CVmix_vars%Vy_cntr(:) = val

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

!BOP

! !IROUTINE: cvmix_att_name
! !INTERFACE:

  function cvmix_att_name(varname)

! !DESCRIPTION:
!  Given a variable short name, returns the precise name of the desired
!  attribute in the cvmix\_data\_type structure.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname

! !OUTPUT PARAMETERS:
    character(len=cvmix_strlen) :: cvmix_att_name

!EOP
!BOC

    select case(trim(varname))
      ! Scalars
      case ("nlev", "NumberLevels", "NumberOfLevels")
        cvmix_att_name = "nlev"
      case ("depth", "ocn_depth", "OceanDepth", "DepthOfOcean")
        cvmix_att_name = "OceanDepth"
      case ('BoundaryLayerDepth','OBL_depth')
        cvmix_att_name = "BoundaryLayerDepth"
      case ("SSH", "surf_hgt", "SeaSurfaceHeight", "SurfaceHeight", "height")
        cvmix_att_name = "SeaSurfaceHeight"
      case ("surf_fric", "SurfaceFriction")
        cvmix_att_name = "SurfaceFriction"
      case ("surf_buoy", "SurfaceBuoyancy", "SurfaceBuoyancyForcing")
        cvmix_att_name = "SurfaceBuoyancyForcing"
      case ("lat", "latitude", "Latitude")
        cvmix_att_name = "lat"
      case ("lon", "longitude", "Longitude")
        cvmix_att_name = "lon"
      case ("coriolis", "Coriolis", "CoriolisFreq", "CoriolisFrequency")
        cvmix_att_name = "Coriolis"
      case ("kOBL_depth", "BoundaryLayerDepthIndex")
        cvmix_att_name = "kOBL_depth"

      ! Variables on level interfaces
      case ("zw", "zw_iface")
        cvmix_att_name = "zw_iface"
      case ("dzw", "dzw_iface")
        cvmix_att_name = "dzw"
      case ("Mdiff", "Udiff", "MomentumDiff", "MomentumDiffusivity")
        cvmix_att_name = "Mdiff_iface"
      case ("Tdiff", "TempDiff", "TemperatureDiff", "TemperatureDiffusivity")
        cvmix_att_name = "Tdiff_iface"
      case ("Sdiff", "SaltDiff", "SalinityDiff", "SalinityDiffusivity")
        cvmix_att_name = "Sdiff_iface"
      case ("Ri", "Ri_iface", "Richardson", "ShearRichardson",                &
            "RichardsonNumber", "ShearRichardsonNumber",                      &
            "ShearRichardson_iface")
        cvmix_att_name = "ShearRichardson_iface"
      case ("buoy", "buoy_iface", "N", "Nsqr", "BuoyancyFreq", "SqrBuoyancy", &
            "SqrBuoyancyFreq", "SqrBuoyancyFreq_iface")
        cvmix_att_name = "SqrBuoyancyFreq_iface"
      case ("kpp_transport", "kpp_nonlocal", "nonlocal_transport",            &
            "nonlocal", "kpp_nonlocal_iface")
        ! Note: this isn't an attribute in the data type, but put / get
        !       uses this as short hand for "both Tnonlocal and Snonlocal"
        cvmix_att_name = "kpp_nonlocal_iface"
      case ("Tnonlocal", "KPP_T_Nonlocal", "kpp_Tnonlocal", "kpp_Ttransport", &
            "kpp_Tnonlocal_iface")
        cvmix_att_name = "kpp_Tnonlocal_iface"
      case ("Snonlocal", "KPP_S_Nonlocal", "kpp_Snonlocal", "kpp_Stransport", &
            "kpp_Snonlocal_iface")
        cvmix_att_name = "kpp_Snonlocal_iface"

      ! Variables on level centers
      case ("z","zt","zt_cntr")
        cvmix_att_name = "zt_cntr"
      case ("dz", "dzt", "CellThickness")
        cvmix_att_name = "dzt"
      case ("rho", "dens", "WaterDensity", "WaterDensity_cntr")
        cvmix_att_name = "WaterDensity_cntr"
      case ("rho_lwr", "dens_lwr", "AdiabWaterDensity",                       &
            "AdiabWaterDensity_cntr")
        cvmix_att_name = "AdiabWaterDensity_cntr"
      case ("Rib", "Ri_bulk", "BulkRichardson", "BulkRichardsonNumber",       &
            "BulkRichardson_cntr")
        cvmix_att_name = "BulkRichardson_cntr"
      case ("Rrho", "strat_param")
        ! Note: this isn't an attribute in the data type, but the I/O routines
        !       use it to denote strat_param_num / strat_param_denom
        cvmix_att_name = "strat_param"
      case ("Rrho_num", "strat_param_num")
        cvmix_att_name = "strat_param_num"
      case ("Rrho_denom", "strat_param_denom")
        cvmix_att_name = "strat_param_denom"
      case ("Buoyancy","buoyancy","buoyancy_cntr")
        cvmix_att_name = "buoyancy_cntr"
      case ("U", "Vx", "Vx_cntr")
        cvmix_att_name = "Vx_cntr"
      case ("V", "Vy", "Vy_cntr")
        cvmix_att_name = "Vy_cntr"
      case DEFAULT
        print*, "ERROR: ", trim(varname), " is not tied to an attribute of",  &
                "the cvmix_data_type structure."
        stop 1
    end select

!EOC

  end function cvmix_att_name

end module cvmix_put_get

