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

  ; Check for memcopy output first
  file_names = "data_"+(/"memcopy", "pointer"/)
  filefound = False
  count = 0
  do while ((.not.filefound).and.(count.lt.2))

    if (check_for_file(file_names(count)+".nc")) then
      print((/"Reading netCDF data from "+file_names(count)/))
      f = addfile(file_names(count)+".nc", "r")
      filefound = True
      x = transpose(f->Tdiff)
      y = f->zw
    else
      if (check_for_file(file_names(count)+".out")) then
        print((/"Reading ascii data from "+file_names(count)/))
        nml  = asciiread("input.nl", -1, "integer")
        nlev = nml(0)
        x = new((/2, nlev+1/), "double")
        y = new(nlev+1, "double")

        data = asciiread(file_names(count)+".out", (/nlev+1,3/), "double")
        filefound = True
        x(0,:) = data(:,1)
        x(1,:) = data(:,2)
        y = data(:,0)
      end if
    end if
    count = count+1
  end do
  if (.not.filefound)
    print((/"ERROR: can not find output; looking for any of following:"/))
    print((/file_names+".nc"/))
    print((/file_names+".out"/))
    exit
  end if

  wks = gsn_open_wks(out_type, "Bryan-Lewis_CVMix")
  gsn_define_colormap(wks, my_color_map)

  ; Basic Graphics set up (don't draw / advance frame to add legend!)
  res = CVMix_res(0, 1.4, -5500, 0)

  ; line & marker styles / axes ranges (y decreases to bottom)
  res@xyMonoMarkLineMode = True
  res@xyMarkLineMode     = "MarkLines"
  res@xyMarkers          = (/0,4/)
  res@xyDashPatterns     = (/0,0/)

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
