load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

; Create Color Table
  my_color_map = (/"White", "Black"/)

  if (isfilepresent("data.nc")) then
    print((/"Reading netCDF data"/))
    f = addfile("data.nc", "r")
    x = f->depth
    y = f->visc
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
  res@gsnDraw  = False
  res@gsnFrame = False

  ; line & marker styles / axes ranges (y decreases to bottom)
  res@xyMonoMarkLineMode = True
  res@xyMarkLineMode     = "Lines"
  res@xyDashPattern      = 0
  res@trXMinF            = 0
  res@trXMaxF            = 1
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
  res@tmXBMajorOutwardLengthF = 0
  res@tmYLMajorOutwardLengthF = 0

  ; Title / axes labels
  res@tiMainString  = "Shear instability"
  res@tiXAxisString = "Local gradient Richardson number"
  res@tiYAxisString = "Normalized viscosity"

  plot = gsn_csm_xy(wks, x, y, res)

  draw(plot)
  frame(wks)

end
