load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

; This script reads in output from the CVMix stand-alone driver using the shear
; mixing technique from Large, et al., 1992. This script takes the output from
; CVMix and uses it to recreate Figure 3 from that paper (page 373).

begin

  ; Create Color Table
  my_color_map = (/"White", "Black"/)

  if (isfilepresent("data.nc")) then
    print((/"Reading netCDF data"/))
    f = addfile("data.nc", "r")
    x = f->Ri
    y = f->diff
  else
    if (isfilepresent("data.out")) then
      print((/"Reading ascii data"/))
      nml  = asciiread("input.nl", -1, "integer")
      nlev = nml(0)
      x = new(nlev+1, "double")
      y = new(nlev+1, "double")

      data = asciiread("data.out", (/nlev+1,2/), "double")
      x = data(:,0)
      y = data(:,1)
      else
        print((/"ERROR: can not find output (looking for data.out or data.nc)"/))
        exit
    end if
  end if

  out_type = "ps"
  ; out_type = "X11"
  wks = gsn_open_wks(out_type, "shear_CVMix")
  gsn_define_colormap(wks, my_color_map)

  ; Basic Graphics set up (don't draw / advance frame to add legend!)
  res          = True
  res@tiMainFuncCode = ":"
  res@tiXAxisFuncCode = res@tiMainFuncCode
  res@tiYAxisFuncCode = res@tiMainFuncCode
  res@gsnDraw  = False
  res@gsnFrame = False

  ; line & marker styles / axes ranges (y decreases to bottom)
  res@xyMonoMarkLineMode = True
  res@xyMarkLineMode     = "Lines"
  res@xyLineThicknessF   = 2.
  res@xyDashPattern      = 0
  res@trXMinF            = 0
  res@trXMaxF            = 0.8
  res@trYMinF            = 0
  res@trYMaxF            = 1

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
  res@tmXBMode = "Manual"
  res@tmXBTickStartF = 0.0
  res@tmXBTickEndF = 0.8
  res@tmXBTickSpacingF = 0.1
  res@tmXBLabelStride = 2
  res@tmYLMode = "Manual"
  res@tmYLTickStartF = 0.0
  res@tmYLTickEndF = 1.0
  res@tmYLTickSpacingF = 0.1
  res@tmYLLabelStride = 2

  ; Title / axes labels
  res@tiMainString  = "Shear Instability"
  res@tiXAxisString = "Ri:B:g:E:"
  res@tiYAxisString = ":F8:k/k:F10::B:0:E:"

  plot = gsn_csm_xy(wks, x, y/5e-3, res)

  draw(plot)
  frame(wks)

end
