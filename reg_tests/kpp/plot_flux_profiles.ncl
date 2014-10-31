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

  if (check_for_file("test3.nc")) then
    print((/"Reading netCDF data"/))
    f = addfile("test3.nc", "r")
    sigma = f->data(0,:)
    phi_m   = f->data(1,:)
    phi_s   = f->data(2,:)
  else
    if (check_for_file("test3.out")) then
      print((/"Reading ascii data"/))
      nlev = 220
      sigma = new((/nlev+1/), "double")
      phi_m   = new((/nlev+1/), "double")
      phi_s   = new((/nlev+1/), "double")

      data = asciiread("test3.out", (/nlev+1,3/), "double")
      sigma = data(:,0)
      phi_m   = data(:,1)
      phi_s   = data(:,2)
      else
        print((/"ERROR: can not find output (looking for test3.out or test3.nc)"/))
        exit
    end if
  end if
  ; We need to know when sigma is positive to stretch the axes correctly and to
  ; be able to plot phi_s only in negative regime.
  ; Use tolerance here because of the way zeta is computed leads to round-off
  ; level values rather than truely 0 values.
  pos_sigma = ind(sigma.gt.-1e-10)
  neg_sigma = ind(sigma.lt.1e-10)

  print((/"For positive sigma values, phi_m and phi_s should be identical."/))
  print((/"Max difference: "+max(abs(phi_m(pos_sigma)-phi_s(pos_sigma)))/))

  ; Compare to LMD94 Figure B1
  wks = gsn_open_wks(out_type, "KPP-flux_profile")
  gsn_define_colormap(wks, my_color_map)

  ; x-axis range will be set before each of two panels is drawn
  res = CVMix_res(0.,0.,0,2)

  ; line & marker styles / axes ranges (y decreases to bottom)
  res@xyMonoMarkLineMode = True
  res@xyMarkLineMode     = "Lines"
  res@xyDashPatterns     = (/0,2/)

  ; Plot / font size, tick marks
  res@tiMainFontHeightF       = 0.02
  res@tiXAxisFontHeightF      = 0.015
  res@tiYAxisFontHeightF      = 0.015
  res@tmXBLabelFontHeightF    = 0.015
  res@tmYLLabelFontHeightF    = 0.015
  res@tmXBMode                = "Explicit"
  res@tmXBMinorOn             = False
  res@tmXBMajorOutwardLengthF = 0
  res@tmYLMode                = "Explicit"
  res@tmYLMinorOn             = False
  res@tmYLValues              = (/0.0, 0.5, 1.0, 1.5, 2.0/)
  res@tmYLMajorOutwardLengthF = 0

  res@tmYRMode                = "Explicit"
  res@tmYRMinorOn             = False
  res@tmYRValues              = (/0.0, 0.5, 1.0, 1.5, 2.0/)
  res@tmYRLabels              = (/"", "", "", "", ""/)
  res@tmYRMajorOutwardLengthF = 0


  res@vpHeightF               = 0.75
  res@vpWidthF                = 0.375

  res1 = res
  ; Title / axes labels
  res1@tiMainOffsetXF = 0.16
  res1@tiXAxisOffsetXF = 0.16
  res1@tiMainString    = "Non-dimensional flux profiles"
  res1@tiXAxisFuncCode = "~"
  res1@tiXAxisString   = "~F33~z~F10~ = d/L = ~F33~s~F10~h/L"
  ; Tick Marks
  res1@tmYLLabels = (/"0", "", "1", "", "2"/)
  res1@trXMinF    = -2.0
  res1@trXMaxF    = 0.
  res1@tmXBValues = (/-2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0./)
  res1@tmXBLabels = (/"-2.", "", "", "", "", "-1.", "", "", "", "", "0"/)
  plot1 = gsn_csm_xy(wks, sigma(neg_sigma), (/phi_m(neg_sigma), phi_s(neg_sigma)/), res1)

  res2 = res
  res2@tmYLLabels = (/"", "", "", "", ""/)
  res2@trXMinF            = 0.
  res2@trXMaxF            = 0.2
  res2@tmXBValues = (/0.0, 0.1, 0.2/)
  res2@tmXBLabels = (/"", "0.1", "0.2"/)
  plot2 = gsn_csm_xy(wks, sigma(pos_sigma), phi_m(pos_sigma), res2)

  txres = True
  txres@txFontHeightF = 0.015
  txres@txFuncCode = "~"
  txtplot1 = gsn_add_text(wks, plot2, "~F33~f~F10~~B~m~N~ = ~F33~f~F10~~B~s~N~", 0.125, 1.4, txres)
  txtplot2 = gsn_add_text(wks, plot1, "~F0~UNSTABLE",       -1.0, 1.0, txres)
  txtplot3 = gsn_add_text(wks, plot2, "~F0~STABLE",          0.1, 1.0, txres)
  txtplot4 = gsn_add_text(wks, plot1, "~F33~f~F10~~B~m~N~", -1.5, 0.465, txres)
  txtplot5 = gsn_add_text(wks, plot1, "~F33~f~F10~~B~s~N~", -0.7, 0.225, txres)
  txtplot6 = gsn_add_text(wks, plot1, "~F33~z~F10~~B~m~N~", -0.2, 1.9, txres)
  txtplot7 = gsn_add_text(wks, plot1, "~F33~z~F10~~B~s~N~", -1.0, 1.9, txres)

  pres = True
;  pres@gsnPanelDebug = True
  pres@gsnPanelXF = (/0.186789, 0.5/)
  gsn_panel(wks, (/plot1, plot2/), (/1,2/), pres)
if (False) then
  temp = gsn_attach_plots(plot1, plot2, False, False)
  draw(plot1)
  frame(wks)
end if

end
