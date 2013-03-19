load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; This script reads in output from the CVMix stand-alone driver using double
; diffusion mixing from Large, et al., 1992. This script takes the output from
; CVMix and uses it to recreate Figures 4a and 4b from that paper (page 374).

begin

; Create Color Table
  my_color_map = (/"White", "Black"/)

  ; Need nlev because we don't use full domain for second plot
  nml  = asciiread("input.nl", -1, "integer")
  nlev = nml(0)
  if (isfilepresent("data.nc")) then
    print((/"Reading netCDF data"/))
    f = addfile("data.nc", "r")
    x1 = f->Rrho
    y1 = f->diff(0:nlev-1,:)
  else
    if (isfilepresent("data.out")) then
      print((/"Reading ascii data"/))
      x1 = new((/nlev,2/), "double")
      y1 = new((/nlev,2/), "double")

      data = asciiread("data.out", (/nlev+1,4/), "double")
      x1 = data(0:nlev-1,0:1)
      y1 = data(0:nlev-1,2:3)
      else
        print((/"ERROR: can not find output (looking for data.out or data.nc)"/))
        exit
    end if
  end if

  out_type = "ps"
  ; out_type = "X11"
  wks = gsn_open_wks(out_type, "Double-Diffusion-Salt")
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
  res@trXMinF            = 1.
  res@trXMaxF            = 2.
  res@trYMinF            = 0.
  res@trYMaxF            = 1.

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
  res@tmXBTickStartF = 1.0
  res@tmXBTickEndF = 2.0
  res@tmXBTickSpacingF = 0.1
  res@tmXBLabelStride = 2
  res@tmYLMode = "Manual"
  res@tmYLTickStartF = 0.0
  res@tmYLTickEndF = 1.0
  res@tmYLTickSpacingF = 0.1
  res@tmYLLabelStride = 2

  ; Title / axes labels
  res@tiMainString  = "Diffusivity in Salt Fingering regime"
  res@tiXAxisString = "R:B::F8:r:E:"
  res@tiYAxisString = ":F8:n:F10::S:d:E::B:s:E::F8: /n:F10::B:f:E:"

  plot = gsn_csm_xy(wks, x1(:,0), 1e3*y1(:,0), res)

  draw(plot)
  frame(wks)

  delete(wks)

  wks = gsn_open_wks(out_type, "Double-Diffusion-Diffusive")
  gsn_define_colormap(wks, my_color_map)

  ; line & marker styles / axes ranges
  res@trXMinF            = 1.
  res@trXMaxF            = 10.
  res@trYMinF            = 0.1
  res@trYMaxF            = 100.
  res@xyYStyle           = "Log"
  res@tmXBTickStartF = 1.0
  res@tmXBTickEndF = 10.0
  res@tmXBTickSpacingF = 1.0
  res@tmXBLabelStride = 2
  res@tmYLMode = "Automatic"
  delete(res@tmYLTickStartF)
  delete(res@tmYLTickEndF)
  delete(res@tmYLTickSpacingF)
  delete(res@tmYLLabelStride)

  ; Plot / font size, tick marks
  res@tmYLMinorOn             = True

  ; Title / axes labels
  res@tiMainString  = "Diffusivity in Diffusive Convection Regime"
  res@tiXAxisString = "1/R:B::F8:r:E:"
  res@tiYAxisString = ":F8:n:F10::S:d:E::F8::B:Q:E: /n"

  plot = gsn_csm_xy(wks, 1.d/x1(0:nlev-1,1), y1(0:nlev-1,1)/1.5e-6, res)

  draw(plot)
  frame(wks)


end
