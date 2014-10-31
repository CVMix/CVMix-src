load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "../common/plot_settings.ncl"
load "../common/version_consistency.ncl"

begin

  ; out_type = "ps"
  ; out_type = "X11"
  out_type = "pdf"

; Create Color Table
  my_color_map = (/"White", "Black"/)

  if (check_for_file("test5.nc")) then
    print((/"Reading netCDF data"/))
    f = addfile("test5.nc", "r")
    zt      = f->zt
    B       = f->buoyancy
    U       = f->U
    Ri_bulk = f->BulkRichardson
    hd      = -f@OBL_depth
  else
    if (check_for_file("test5.out")) then
      print((/"Reading ascii data"/))
      nlev = 10
      zt      = new((/nlev/), "double")
      B       = new((/nlev/), "double")
      U       = new((/nlev/), "double")
      Ri_bulk = new((/nlev/), "double")

      data = asciiread("test5.out", (/nlev,4/), "double")
      zt      = data(:,0)
      Ri_bulk = data(:,1)
      U       = data(:,2)
      B       = data(:,3)
      hd      = new(1,"double")
      else
        print((/"ERROR: can not find output (looking for test5.out or test5.nc)"/))
        exit
    end if
  end if
  true_zt = ispan(-50,0,1)*1.d
  true_B  = new((/51/), "double")
  true_U  = new((/51/), "double")
  true_Ri = new((/51/), "double")
  true_Vt = new((/51/), "double")
  ; Params for true_Vt
  Cv  = 1.5d
  Bt  = -0.2d
  cs  = 24.d*sqrt(17.d)
  eps = 0.1d
  N   = 0.01d
  Ric = 0.3d
  kap = 0.4d
  w_s = 0.01d * kap
  do i=0,50
    if (true_zt(i).lt.-17.d) then
      true_B(i)  = 0.0001d + 0.0001d*(true_zt(i) + 17.d)
      true_U(i)  = 0.1d +0.1d*(true_zt(i) + 17.d)/33.d
      true_Vt(i) = Cv*sqrt(-Bt/(cs*eps))*(-true_zt(i))*N*w_s/(Ric*kap^2)
      true_Ri(i) = -true_zt(i)*(0.0001d-true_B(i)) /                          \
                   ((0.1d-true_U(i))^2 + true_Vt(i))
    else
      true_B(i)  = 0.0001d
      true_U(i)  = 0.1d
      true_Ri(i) = 0.d
    end if
  end do

  ; Want to make three plots, all with depth on y-axis (min: -50, max: 0)
  ; 1: B on x-axis (min: -0.003, max: 0.001)
  ; 2: U on x-axis (min: 0, max: 0.12)
  ; 3: Ri_bulk on x-axis (min: 0, max: 3)

  ; Should be identical to LMD94 Figure C1
  wks = gsn_open_wks(out_type, "KPP-bulk_Rich")
  gsn_define_colormap(wks, my_color_map)

  ; x-axis range will be set before each of three panels is drawn
  res = CVMix_res(0.d,0.d,-50.d,0.d)

  ; no label on x-axis
  res@tiXAxisString = ""
  res@tiYAxisString = ""

  ; Plot / font size, tick marks
  res@tiMainFontHeightF       = 0.02
  res@tiXAxisFontHeightF      = 0.015
  res@tiYAxisFontHeightF      = 0.015
  res@tmXTLabelFontHeightF    = 0.015
  res@tmYLLabelFontHeightF    = 0.015
  res@tmXUseBottom            = False
  res@tmXBMode                = "Explicit"
  res@tmXTMode                = "Explicit"
  res@tmXTLabelsOn            = True
  res@tmXBMinorOn             = False
  res@tmXBMajorOutwardLengthF = 0
  res@tmXTMinorOn             = False
  res@tmXTMajorOutwardLengthF = 0
  res@tmXTLabelFont           = "duplex_roman"
  res@tmYLMode                = "Explicit"
  res@tmYLMinorOn             = False
  res@tmYLMajorOutwardLengthF = 0
  res@tmYLLabelFont           = "duplex_roman"
  res@tmYLValues              = (/-50, -45, -40, -35, -30, -25, -20, -15,     \
                                  -10,  -5,   0/)
  res@vpHeightF               = 0.7

  plot  = new(3,graphic)
  plot2 = new(3,graphic)

  ; Buoyancy plot: common settings
  bres = res
  bres@tiMainString  = "Buoyancy (m/s:S:2:E: )"
  bres@vpWidthF      = 0.2
  bres@trXMinF       = -0.003d
  bres@trXMaxF       = 0.001d
  bres@tmXBValues    = (/-0.003d, -0.002d, -0.001d, 0.d, 0.001d/)
  bres@tmXTValues    = bres@tmXBValues
  bres@tmXBLabels    = (/     "", "", "",  "", ""/)
  bres@tmXTLabels    = (/"-.003", "", "", "0", ""/)
  bres@tiYAxisString = "depth (m)"
  bres@tmYLLabels    = (/"-50", "", "-40", "", "-30", "", "-20", "", "-10", "", "0"/)

  ; Plot 1: Markers for discrete values
  bres@xyMarkLineMode = "Markers"
  bres@xyMarker       = 7 ; triangle
  plot(0) = gsn_csm_xy(wks, B, zt, bres)
  ; Plot 2: Line for continues values
  bres@xyMarkLineMode = "Lines"
  plot2(0) = gsn_csm_xy(wks, true_B, true_zt, bres)
  ; Combine into one plot
  overlay(plot(0), plot2(0))

  ; Velocity plot: common settings
  ures = res
  ures@tiMainString  = "U (m/s)"
  ures@vpWidthF      = 0.2
  ures@trXMinF       = 0.d
  ures@trXMaxF       = 0.12d
  ures@tmXBValues    = (/0.d, 0.02d, 0.04d, 0.06d, 0.08d, 0.1d/) 
  ures@tmXTValues    = ures@tmXBValues
  ures@tmXBLabels    = (/ "", "", "", "",  "",   ""/)
  ures@tmXTLabels    = (/"0", "", "", "", "", "0.1"/)
  ures@tiYAxisString = "" 
  ures@tmYLLabels    = (/"", "", "", "", "", "", "", "", "", "", ""/)

  ; Plot 1: Markers for discrete values
  ures@xyMarkLineMode = "Markers"
  ures@xyMarker       = 7 ; triangle
  plot(1) = gsn_csm_xy(wks, U, zt, ures)
  ; Plot 2: Line for continues values
  ures@xyMarkLineMode = "Lines"
  plot2(1) = gsn_csm_xy(wks, true_U, true_zt, ures)
  ; Combine into one plot
  overlay(plot(1), plot2(1))

  ; bulk Richardson plot: common settings
  rres = res
  rres@tiMainString  = "bulk Richardson"
  rres@vpWidthF      = 0.3
  rres@trXMinF       = 0.d
  rres@trXMaxF       = 3.d
  rres@tmXBValues    = (/0.d, 1.d, 2.d/) 
  rres@tmXTValues    = (/0.d, 0.3d, 1.d, 2.d/) 
  rres@tmXBLabels    = (/  "",              "",  ""/)
  rres@tmXTLabels    = (/ "0", "Ri:B:c:E:", "", "2"/)
  rres@tiYAxisString = ""
  rres@tmYLLabels    = (/"", "", "", "", "", "", "", "", "", "", ""/)

  ; Plot 1: Markers for discrete values
  rres@xyMarkLineMode = "Markers"
  rres@xyMarker       = 7 ; triangle
  plot(2) = gsn_csm_xy(wks, Ri_bulk, zt, rres)
  ; Plot 2: Line for continues values
  rres@xyMarkLineMode = "Lines"
  plot2(2) = gsn_csm_xy(wks, true_Ri, true_zt, rres)
  ; Combine into one plot
  overlay(plot(2), plot2(2))

  ; Add some text
  ; (1) label for h_c
  rres@txFont         = "duplex_roman"
  rres@txFontHeightF  = 0.01
  plot_text = gsn_add_text(wks, plot(2), "- h:B:c:E:", Ri_bulk(3)+1.2d*Ric,   \
                                zt(3), rres)

  ; (2) Marker and label for h_d
  ; If output is netCDF, place an 'x' on the discrete OBL depth
  ; (OBL_depth is not saved in ascii output)
  if (.not.ismissing(hd)) then
    plot_text2 = gsn_add_text(wks, plot(2), "x", Ric, hd, rres)
    plot_text3 = gsn_add_text(wks, plot(2), "- h:B:d:E:", 2.2d*Ric, hd, rres)
  end if

  ; Combine all plots into single page
  pan_res = True
;  pan_res@gsnPanelDebug = True
  gsn_panel(wks,plot,(/1,4/),pan_res)
  ; For some reason setting up a panel with just 3 plots causes the automatic 
  ; resizing to push the third panel off the page, but 4 plots are okay.

end
