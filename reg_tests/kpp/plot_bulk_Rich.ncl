load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

; Create Color Table
  my_color_map = (/"White", "Black"/)

  if (isfilepresent("test5.nc")) then
    print((/"Reading netCDF data"/))
    f = addfile("test5.nc", "r")
    zt      = f->zt
    B       = f->buoyancy
    U       = f->Vx
    Ri_bulk = f->Ri_bulk
  else
    if (isfilepresent("test5.out")) then
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
      else
        print((/"ERROR: can not find output (looking for test5.out or test5.nc)"/))
        exit
    end if
  end if
  true_zt = ispan(-50,0,1)*1.d
  true_B = new((/51/), "double")
  do i=0,50
    if (true_zt(i).lt.-17.d) then
      true_B(i) = 0.0001d + 0.0001d*(true_zt(i) + 17.d)
    else
      true_B(i) = 0.0001d
    end if
  end do

  ; out_type = "ps"
  out_type = "X11"

  ; Want to make three plots, all with depth on y-axis (min: -50, max: 0)
  ; 1: B on x-axis (min: -0.003, max: 0.001)
  ; 2: U on x-axis (min: 0, max: 0.12)
  ; 3: Ri_bulk on x-axis (min: 0, max: 3)

  wks = gsn_open_wks(out_type, "KPP-bulk_Rich")
  gsn_define_colormap(wks, my_color_map)

  ; Basic Graphics set up (don't draw / advance frame to add legend!)
  res          = True
  res@gsnDraw  = False
  res@gsnFrame = False
  res@gsnMaximize = True

  ; line & marker styles / axes ranges (y decreases to bottom)
  res@xyMarkLineMode     = "Markers"
  res@xyMarker           = 7 ; triangle
  res@trXMinF            = -0.003d
  res@trXMaxF            = 0.001d
  res@trYMinF            = -50.d
  res@trYMaxF            = 0.d

  ; Plot / font size, tick marks
;  res@vpHeightF               = 0.75
;  res@vpWidthF                = 0.75
  res@tiMainFontHeightF       = 0.02
  res@tiXAxisFontHeightF      = 0.015
  res@tiYAxisFontHeightF      = 0.015
  res@tmXBLabelFontHeightF    = 0.015
  res@tmYLLabelFontHeightF    = 0.015
  res@tmXBMode                = "Explicit"
  res@tmXTMode                = "Explicit"
  res@tmXBMinorOn             = False
  res@tmXBValues              = (/-0.003d, -0.002d, -0.001d, 0.d, 0.001d/)
  res@tmXTValues              = res@tmXBValues
  res@tmXBLabels              = (/     "", "", "",  "", ""/)
  res@tmXTLabels              = (/"-.003", "", "", "0", ""/)
  res@tmXBMajorOutwardLengthF = 0
  res@tmYLMode                = "Explicit"
  res@tmYLMinorOn             = False
  res@tmYLValues              = (/-50, -45, -40, -35, -30, -25, -20, -15,     \
                                  -10,  -5,   0/)
  res@tmYLLabels              = (/"50", "", "40", "", "30", "", "20", "",     \
                                  "10", "",  "0"/)
  res@tmYLMajorOutwardLengthF = 0

  ; Title / axes labels
  res@tiMainString  = "Buoyancy (m/s:S:2:E: )"
  res@tiXAxisString = ""
  res@tiYAxisString = "depth (m)"

  plot = gsn_csm_xy(wks, B, zt, res)

  res@xyMarkLineMode     = "Lines"
  plot2 = gsn_csm_xy(wks, true_B, true_zt, res)

  overlay(plot, plot2)
  draw(plot)
  frame(wks)

end
