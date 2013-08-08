load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

; Create Color Table
  my_color_map = (/"White", "Black"/)

  if (isfilepresent("test3.nc")) then
    print((/"Reading netCDF data"/))
    f = addfile("test3.nc", "r")
    sigma = f->data(0,:)
    w_m   = f->data(1,:)
    w_s   = f->data(2,:)
  else
    if (isfilepresent("test3.out")) then
      print((/"Reading ascii data"/))
      nlev = 220
      sigma = new((/nlev+1/), "double")
      w_m   = new((/nlev+1/), "double")
      w_s   = new((/nlev+1/), "double")

      data = asciiread("test3.out", (/nlev+1,3/), "double")
      sigma = data(:,0)
      w_m   = data(:,1)
      w_s   = data(:,2)
      else
        print((/"ERROR: can not find output (looking for test3.out or test3.nc)"/))
        exit
    end if
  end if
  pos_sigma = ind(sigma.gt.0)
  sigma(pos_sigma) = 10.d*sigma(pos_sigma)

  ; out_type = "ps"
  out_type = "X11"
  wks = gsn_open_wks(out_type, "KPP-flux_profile")
  gsn_define_colormap(wks, my_color_map)

  ; Basic Graphics set up (don't draw / advance frame to add legend!)
  res          = True
  res@gsnMaximize = True

  ; line & marker styles / axes ranges (y decreases to bottom)
  res@xyMonoMarkLineMode = True
  res@xyMarkLineMode     = "Lines"
  res@xyDashPatterns     = (/0,2/)
  res@trXMinF            = -2.0
  res@trXMaxF            = 2.0
  res@trYMinF            = 0
  res@trYMaxF            = 2

  ; Plot / font size, tick marks
;  res@vpHeightF               = 0.75
;  res@vpWidthF                = 0.75
  res@tiMainFontHeightF       = 0.02
  res@tiXAxisFontHeightF      = 0.015
  res@tiYAxisFontHeightF      = 0.015
  res@tmXBLabelFontHeightF    = 0.015
  res@tmYLLabelFontHeightF    = 0.015
  res@tmXBMode                = "Explicit"
  res@tmXBMinorOn             = False
  res@tmXBValues              = (/-2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8,   \
                                  -0.6, -0.4, -0.2,  0.0,  1.0,  2.0/)
  res@tmXBLabels              = (/"-2.", "", "", "", "", "-1.", "", "", "",   \
                                  "", "0", ".1", ".2"/)
  res@tmXBMajorOutwardLengthF = 0
  res@tmYLMode                = "Explicit"
  res@tmYLMinorOn             = False
  res@tmYLValues              = (/0.0, 0.5, 1.0, 1.5, 2.0/)
  res@tmYLLabels              = (/"0", "", "1", "", "2"/)
  res@tmYLMajorOutwardLengthF = 0

  ; Title / axes labels
  res@tiMainString  = ""
  res@tiXAxisString = "Zeta = d/L = sigma h/L"
  res@tiYAxisString = ""

  plot = gsn_csm_xy(wks, sigma, (/w_m, w_s/), res)

end
