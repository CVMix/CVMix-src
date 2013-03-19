load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

; Create Color Table
  my_color_map = (/"White", "Black"/)

  if (isfilepresent("data.nc")) then
    print((/"Reading netCDF data"/))
    f = addfile("data.nc", "r")
    x1 = f->Rrho(:,0)
    y1 = f->diff(:,0)
  else
    if (isfilepresent("data.out")) then
      print((/"Reading ascii data"/))
      nml  = asciiread("input.nl", -1, "integer")
      nlev = nml(0)
      x1 = new((/nlev+1,2/), "double")
      y1 = new((/nlev+1,2/), "double")

      data = asciiread("data.out", (/nlev+1,4/), "double")
      x1 = data(:,0:1)
      y1 = data(:,2:3)
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
  res@gsnDraw  = False
  res@gsnFrame = False

  ; line & marker styles / axes ranges (y decreases to bottom)
  res@xyMonoMarkLineMode = True
  res@xyMarkLineMode     = "Lines"
  res@xyDashPattern      = 0
  res@trXMinF            = 1
  res@trXMaxF            = 2
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
  res@tiMainString  = "Diffusivity in Salt Fingering regime"
  res@tiXAxisString = "Rrho"
  res@tiYAxisString = "Normalized diffusivity"

  plot = gsn_csm_xy(wks, x1(:,0), 1e3*y1(:,0), res)

  draw(plot)
  frame(wks)

  delete(res)
  delete(wks)

  wks = gsn_open_wks(out_type, "Double-Diffusion-Diffusive")
  gsn_define_colormap(wks, my_color_map)

  ; Basic Graphics set up (don't draw / advance frame to add legend!)
  res          = True
  res@gsnDraw  = False
  res@gsnFrame = False

  ; line & marker styles / axes ranges (y decreases to bottom)
  res@xyMonoMarkLineMode = True
  res@xyMarkLineMode     = "Lines"
  res@xyDashPattern      = 0
  res@trXMinF            = 1
  res@trXMaxF            = 10
  res@trYMinF            = 0.1
  res@trYMaxF            = 100
  res@xyYStyle           = "Log"

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
  res@tiMainString  = "Diffusivity in Diffusive Convection Regime"
  res@tiXAxisString = "1 / Rrho"
  res@tiYAxisString = "Normalized diffusivity"

  plot = gsn_csm_xy(wks, 1.d/x1(0:nlev-1,1), y1(0:nlev-1,1)/1.5e-6, res)

  draw(plot)
  frame(wks)


end
