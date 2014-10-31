load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "../common/plot_settings.ncl"
load "../common/version_consistency.ncl"

begin

  out_type = "pdf"
  ; out_type = "ps"
  ; out_type = "X11"

; Create Color Table
  my_color_map = (/"White", "Black"/)

  in_file = "single_col.nc"
  if (check_for_file(in_file)) then
    print((/"Reading netCDF data"/))
    f = addfile(in_file, "r")
    x = f->Tdiff
    y = f->zw
    nlev = dimsizes(y)
  else
    print((/"ERROR: can not find output (looking for "+in_file+")"/))
    exit
  end if

  wks = gsn_open_wks(out_type, "Simmons-single_col")
  gsn_define_colormap(wks, my_color_map)

  ; Basic Graphics set up (don't draw / advance frame to add legend!)
  res = CVMix_res(0,1e4*max(x)*1.05,min(y),0)
  res@gsnDraw  = False
  res@gsnFrame = False

  ; line & marker styles / axes ranges (y decreases to bottom)
  res@xyLineThicknessF   = 2

  ; Plot / font size, tick marks
  res@vpHeightF               = 0.55
  res@vpWidthF                = 0.77
  res@tiMainFontHeightF       = 0.02
  res@tiXAxisFontHeightF      = 0.015
  res@tiYAxisFontHeightF      = 0.015
  res@tmXBLabelFontHeightF    = 0.015
  res@tmYLLabelFontHeightF    = 0.015
  res@tmXBMinorOn             = False
  res@tmYLMinorOn             = False
  res@tmXBMajorOutwardLengthF = 0
  res@tmYLMajorOutwardLengthF = 0

  ; Title / axes labels
  res@tiMainString  = "Vertical tracer diffusivity due to tidal energy"
  res@gsnCenterString = "Location is ("+f@column_lon+", "+f@column_lat+")"
  res@tiXAxisString = "Diffusivity (1e-4 m:S:2:N:/sec)"
  res@tiYAxisString = "Depth (m)"

  print((/"Making plot..."/))
  print((/res@gsnCenterString/))
  plot = gsn_csm_xy(wks, 1e4*x(0:nlev-2), y(0:nlev-2), res)

  draw(plot)
  frame(wks)

end
