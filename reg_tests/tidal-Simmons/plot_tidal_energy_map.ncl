load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

  in_file = "../../inputdata/tidal_energy_gx1v6_20130512.nc"
  if (isfilepresent(in_file)) then
    print((/"Reading netCDF data"/))
    f = addfile(in_file, "r")
    tmplat = f->lat
    tmplon = f->lon
    dims = dimsizes(tmplat)
    nlat = dims(0)
    nlon = dims(1)+1
    lon = new((/nlat, nlon/), typeof(tmplon))
    lat = new((/nlat, nlon/), typeof(tmplat))
    lat(:,1:nlon-1) = tmplat
    lon(:,1:nlon-1) = tmplon
    lat(:,0) = tmplat(:,nlon-2)
    lon(:,0) = tmplon(:,nlon-2)

    ; data = log10(f->TIDAL_ENERGY_FLUX(1:nlat-1,:))
    data = f->TIDAL_ENERGY_FLUX(1:nlat-1,:)
    data@lat2d = lat
    data@lon2d = lon
  else
    print((/"ERROR: can not find output (looking for "+in_file+")"/))
    exit
  end if

  out_type = "ps"
  ; out_type = "X11"
  wks = gsn_open_wks(out_type, "tidal_energy_flux")
  gsn_define_colormap(wks, "BlAqGrYeOrRe")

  ; Basic Graphics set up (don't draw / advance frame to add legend!)
  res          = True
;  res@mpCenterLonF = -180.
  res@cnFillOn = True
  res@cnFillMode = "RasterFill"
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels = 10^(0.5*ispan(-9,1,1))
  res@lbLabelAlignment = "ExternalEdges"
  res@lbLabelStrings = (/"10:S:-5:E:","","10:S:-4:E:","","10:S:-3:E:","","10:S:-2:E:","","10:S:-1:E:","","1","","10"/)
  res@cnLinesOn = False
  res@gsnSpreadColors = True
;  res@gsnSpreadColorStart = 10
;  res@gsnSpreadColorEnd   = 96

  ; Title / axes labels
  res@tiMainString  = "Tidal Energy Flux"
;  res@tiXAxisString = "Diffusivity (1e-4 m:S:2:N:/sec)"
;  res@tiYAxisString = "Depth (m)"

  print((/"Making Plot"/))
  plot = gsn_csm_contour_map_ce(wks, data, res)

end
