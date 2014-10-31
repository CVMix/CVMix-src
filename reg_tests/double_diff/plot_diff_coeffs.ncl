load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "../common/plot_settings.ncl"
load "../common/version_consistency.ncl"

; This script reads in output from the CVMix stand-alone driver using double
; diffusion mixing from Large, et al., 1992. This script takes the output from
; CVMix and uses it to recreate Figures 4a and 4b from that paper (page 374).

begin

  out_type = "pdf"
  ; out_type = "ps"
  ; out_type = "X11"

; Create Color Table
  my_color_map = (/"White", "Black"/)

  ; Need nlev because we don't use full domain for second plot
  nml  = asciiread("input.nl", -1, "integer")
  nlev = 2*nml(0)
  yt   = new((/2, nlev/), "double")
  if (check_for_file("data.nc")) then
    print((/"Reading netCDF data"/))
    f = addfile("data.nc", "r")

    x1      = f->Rrho
    ys      = f->Sdiff(0:nlev-1,0)
    yt(0,:) = f->Tdiff(0:nlev-1,0)
    yt(1,:) = f->Tdiff(0:nlev-1,1)
  else
    if (check_for_file("data.out")) then
      print((/"Reading ascii data"/))
      x1 = new((/nlev  /), "double")
      ys = new((/nlev  /), "double")

      data = asciiread("data.out", (/nlev+1,5/), "double")
      x1      = data(0:nlev-1,  0)
      ys      = data(0:nlev-1,  3)
      yt(0,:) = data(0:nlev-1,  1)
      yt(1,:) = data(0:nlev-1,  2)
      else
        print((/"ERROR: can not find output (looking for data.out or data.nc)"/))
        exit
    end if
  end if

  wks = gsn_open_wks(out_type, "ddiff-salt")
  gsn_define_colormap(wks, my_color_map)

  ; Basic Graphics set up (don't draw / advance frame to add legend!)
  ; axes ranges are for convection regime
  res = CVMix_res(1., 2., 0., 1.)

  ; line & marker styles for convection regime
  res@xyMonoMarkLineMode = True
  res@xyMarkLineMode     = "Lines"
  res@xyLineThicknessF   = 2.
  res@xyDashPatterns     = (/0,1/)

  ; Plot / font size, tick marks
  res@vpHeightF            = 0.55
  res@vpWidthF             = 0.77
  res@tiMainFontHeightF    = 0.02
  res@tiXAxisFontHeightF   = 0.015
  res@tiYAxisFontHeightF   = 0.015
  res@tmXBLabelFontHeightF = 0.015
  res@tmYLLabelFontHeightF = 0.015
  res@tmXBMinorOn          = False
  res@tmYLMinorOn          = False
  res@tmXBMode             = "Explicit"
  res@tmXBValues           = 0.1*ispan(10,20,1)
  res@tmXBLabels           = (/"1.00", "", "1.20", "", "1.40", "", "1.60", \
                           "", "1.80", ":F10:R:F8::B:r:E::F10::S:0:E:", "2.00"/) 
  res@tmYLMode             = "Explicit"
  res@tmYLValues           = 0.1*ispan(0,10,1)
  res@tmYLLabels           = (/"0.00", "", "0.20", "", "0.40", "",            \
                               "0.60", "", "0.80", "", "1.00"/) 

  ; Title / axes labels
  res@tiMainString  = "Diffusivity in Salt Fingering regime"
  res@tiXAxisString = ":F10:R:B::F8:r:E:"
  res@tiYAxisString = ":F8:k/k:F10::B:0:E:"

  plot = gsn_csm_xy(wks, x1(0:nlev/2-1), ys(0:nlev/2-1)/ys(0), res)

  draw(plot)
  frame(wks)

  delete(wks)

  wks = gsn_open_wks(out_type, "ddiff-diffuse")
  gsn_define_colormap(wks, my_color_map)

  ; line & marker styles / axes ranges for salt fingering regime
  res = CVMix_res(1., 10., 0.1, 100.)

  res@xyYStyle         = "Log"
  delete(res@tmXBValues)
  delete(res@tmXBLabels)
  res@tmXBValues       = ispan(1,10,1)
  res@tmXBLabels       = (/"1.0", "", "3.0", "", "5.0", "", "7.0", "", "9.0", ""/)
  res@tmYLMode         = "Automatic"
  res@tmYLMinorPerMajor = 8
  res@tmYRMinorPerMajor = 8

  ; Plot / font size, tick marks
  res@tmYLMinorOn             = True

  ; Title / axes labels
  res@tiMainString  = "Diffusivity in Diffusive Convection Regime"
  res@tiXAxisString = "1/:F10:R:B::F8:r:E:"
  res@tiYAxisString = ":F8:k/n:F10::B:mol:E:"

  plot = gsn_csm_xy(wks, 1.d/x1(nlev/2:nlev-1), yt(:,nlev/2:nlev-1)/1.5d-6, res)

  draw(plot)
  frame(wks)


end
