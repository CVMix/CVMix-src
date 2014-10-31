load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "../common/plot_settings.ncl"
load "../common/version_consistency.ncl"

; This script reads in output from the CVMix stand-alone driver using the shear
; mixing technique from Large, et al., 1992. This script takes the output from
; CVMix and uses it to recreate Figure 3 from that paper (page 373).

begin

  out_type = "pdf"
  ; out_type = "ps"
  ; out_type = "X11"

  ; Create Color Table
  my_color_map = (/"White", "Black"/)

  ; Read data for PP81 plot
  if (check_for_file("data_PP2d.nc")) then
    print((/"Reading netCDF data"/))
    f = addfile("data_PP2d.nc", "r")
    x = f->ShearRichardson
    y = transpose(f->Mdiff)
  else
    if (check_for_file("data_PP2d.out")) then
      print((/"Reading ascii data"/))
      nml  = asciiread("input.nl", -1, "integer")
      nlev = nml(0)
      x = new(nlev+1, "double")
      y = new((/6,nlev+1/), "double")

      data = asciiread("data_PP2d.out", (/nlev+1,7/), "double")
      x = data(:,0)
      do i=1,6
        y(i-1,:) = data(:,i)
      end do
      else
        print((/"ERROR: can not find output (looking for data_pp2d.out or .nc)"/))
        exit
    end if
  end if

  wks = gsn_open_wks(out_type, "PP_shear-instability")
  gsn_define_colormap(wks, my_color_map)

  ; Basic Graphics set up (don't draw / advance frame to add legend!)
  res = CVMix_res(0,1,0,100)

  ; line & marker styles / axes ranges (y decreases to bottom)
  B = 3
  C = NhlNewDashPattern(wks,"$$$$$$$$___$___$___$___$___$$$$$$$$")
  D = 0
  E = NhlNewDashPattern(wks,"$$$$$$$$__$__$__$$$$$$$$")
  F = 2
  G = NhlNewDashPattern(wks,"$$$$$$__$$$$$$__$$$$$$__$$$$$$__$$$$$$__")
  res@xyMonoMarkLineMode = True
  res@xyMarkLineMode     = "Lines"
  res@xyDashPatterns     = (/B,C,D,E,F,G/)
  res@xyLineThicknessF   = 2.

  ; Plot / font size, tick marks
  res@vpHeightF            = 0.55
  res@vpWidthF             = 0.55
  res@tiMainFontHeightF    = 0.02
  res@tiXAxisFontHeightF   = 0.015
  res@tiYAxisFontHeightF   = 0.015
  res@tmXBLabelFontHeightF = 0.015
  res@tmYLLabelFontHeightF = 0.015
  res@tmXBMinorOn          = False
  res@tmYLMinorOn          = False
  res@tmXBMode             = "Explicit"
  res@tmXBValues           = 0.1*ispan(1,10,1)
  res@tmXBLabels           = (/".1", ".2", ".3", ".4", ".5", ".6", ".7",      \
                               ".8", ".9", "1.0"/)
  res@tmYLMode = "Manual"
  res@tmYLTickStartF = 0.0
  res@tmYLTickEndF   = 100.
  res@tmYLTickSpacingF = 10.
  res@tmYLLabelStride = 1

  ; Title / axes labels
  res@tiMainString  = "Dependence of :F8:n:F25: on Ri"
  res@tiXAxisString = ":F25:Ri"
  res@tiYAxisString = ":F8:n:F25: (cm:S:2:N:/s)"

  plot = gsn_csm_xy(wks, x, y*10000.d, res)

  ; Add text labels
  txres = CVMix_txres(0.015)

  label1 = gsn_add_text(wks, plot, "B", 0.09, 13, txres)
  label1 = gsn_add_text(wks, plot, "C", 0.08, 31, txres)
  label1 = gsn_add_text(wks, plot, "D", 0.30, 19, txres)
  label1 = gsn_add_text(wks, plot, "E", 0.58, 29, txres)
  label1 = gsn_add_text(wks, plot, "F", 0.37, 22, txres)
  label1 = gsn_add_text(wks, plot, "G", 0.22, 18, txres)

  draw(plot)
  frame(wks)

end
