load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

; Create Color Table
  my_color_map = (/"White", "Black"/)

  nml   = asciiread("../exe/low_lat.nl", -1, "integer")
  nlev  = nml(0)
  datal = asciiread("../exe/diff1.txt", (/nlev+1,2/), "double")
  datah = asciiread("../exe/diff2.txt", (/nlev+1,2/), "double")

  x = new((/2, nlev+1/), "double")
  x(0,:) = datal(:,1)
  x(1,:) = datah(:,1)

  y = new(nlev+1, "double")
  y = datal(:,0)

  ; out_type = "ps"
  out_type = "X11"
  wks = gsn_open_wks(out_type, "Bryan-Lewis_CVMix")
  gsn_define_colormap(wks, my_color_map)

  ; Basic Graphics set up (don't draw / advance frame to add legend!)
  res          = True
  res@gsnDraw  = False
  res@gsnFrame = False

  ; line & marker styles / axes ranges (y decreases to bottom)
  res@xyMonoMarkLineMode = True
  res@xyMarkLineMode     = "MarkLines"
  res@xyMarkers          = (/0,4/)
  res@xyDashPatterns     = (/0,0/)
  res@trXMinF            = 0
  res@trXMaxF            = 1.4
  res@trYMinF            = 0
  res@trYMaxF            = 5500
  res@trYReverse         = True

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
  res@tiMainString  = "Background vertical tracer diffusivity"
  res@tiXAxisString = "Diffusivity (1e-4 m:S:2:N:/sec)"
  res@tiYAxisString = "Depth (m)"

  plot = gsn_csm_xy(wks, 1e4*x, y, res)

  ; Legend
  lres = True
  lres@vpWidthF         = 0.30
  lres@vpHeightF        = 0.14
  lres@lgMonoItemType   = True
  lres@lgItemType         = "Markers"
  lres@lgAutoManage     = False
  lres@lgLabelFontHeightF = 0.035
  lres@lgMonoMarkerIndex = False
  lres@lgMarkerIndexes  = (/0, 4/)
  lres@lgLineLabelsOn   = False

  leg_text = (/"Low", "High"/)+" latitudes"
  legend = gsn_create_legend(wks, 2, leg_text, lres)

  ares = True
  ares@amParallelPosF = 0.27
  ares@amOrthogonalPosF = -0.33
  annotate = gsn_add_annotation(plot, legend, ares)
  draw(plot)
  frame(wks)

end
