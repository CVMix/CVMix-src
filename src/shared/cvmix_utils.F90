module cvmix_utils

!BOP
!\newpage
! !MODULE: cvmix_utils
!
! !DESCRIPTION:
!  This module contains routines that are called by multiple modules but don't
!  specifically compute anything mixing related.
!\\
!\\

! !USES:

   use cvmix_kinds_and_types, only : cvmix_r8,                  &
                                     cvmix_strlen
!EOP

  implicit none
  private
  save

!BOP

! !PUBLIC MEMBER FUNCTIONS:
  public :: cvmix_att_name

!EOP

contains

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

end module cvmix_utils
