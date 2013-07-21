load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

  ; out_type = "ps"
  out_type = "pdf"
  ; out_type = "x11"

  N = 2000
  Nvert = 40
  tol = 1e-4
  plot_at = -22.05
  make_last_plot = False
  all_pos = True
  ; all_pos = False
  hmix   = new(N+1, float)
  cont   = new(N+1, float)
  linear = new(N+1, float)
  quad   = new(N+1, float)
  cubic  = new(N+1, float)

  depth = (/0., -10., -20., -30., -40/)
  midpoints = (/-5., -15., -25., -35./)
  Ncells = dimsizes(midpoints)
  v_depth   = new(Nvert*(Ncells-1)+1, float)
  do ii=0,Ncells-2
    slope = (midpoints(ii+1)-midpoints(ii))/Nvert
    do jj=0, Nvert
      kk = ii*Nvert+jj
      v_depth(kk) = slope*jj + midpoints(ii)
    end do
  end do
  Ri_crit   = new(Nvert*(Ncells-1)+1, float)
  Ri_crit   = 0.3
  exact     = new(Nvert*(Ncells-1)+1, float)
  lin_col   = new(Nvert*(Ncells-1)+1, float)
  quad_col  = new(Nvert*(Ncells-1)+1, float)
  cubic_col = new(Nvert*(Ncells-1)+1, float)
  Ri_bulk = new(Ncells+1, float)
  Rib_mid = new(Ncells, float)
  total_its = 0

  do i=0,N
    hmix(i) = -(10. + i*20./N)
    if (abs(hmix(i)-plot_at).lt.1e-6) then
      extra_plot = True
      make_last_plot = True
      icol = i
    else
      extra_plot = False
    end if
    cont(i) = hmix(i)-2.
    Ri_bulk(0)  = 0.
    Ri_bulk(1)  = 0.
    Rib_mid(0) = 0.
    if (hmix(i).gt.-20.) then
      Ri_bulk(2) = 0.15*(20.+hmix(i))
      Ri_bulk(3) = 1.5+Ri_bulk(2)
      Ri_bulk(4) = 1.5+Ri_bulk(3)
      Rib_mid(1) = 0.0075*(20.+hmix(i))*(20.+hmix(i))
      Rib_mid(2) = 0.5*(Ri_bulk(2) + Ri_bulk(3))
      Rib_mid(3) = 0.5*(Ri_bulk(3) + Ri_bulk(4))
    else 
      Ri_bulk(2)  = 0.
      Rib_mid(1) = 0.
      if (hmix(i).gt.-30.) then
        Ri_bulk(3) = 0.15*(30.+hmix(i))
        Ri_bulk(4) = 1.5+Ri_bulk(3)
        Rib_mid(2) = 0.0075*(30.+hmix(i))*(30.+hmix(i))
        Rib_mid(3) = 0.5*(Ri_bulk(3) + Ri_bulk(4))
      else
        Ri_bulk(3) = 0.
        Ri_bulk(4) = 0.15*(40.+hmix(i))
        Rib_mid(2) = 0.
        Rib_mid(3) = 0.0075*(40.+hmix(i))*(40.+hmix(i))
      end if
    end if
    if (Rib_mid(1).gt.0.3) then
      ; Rib(-5) < 0.3, Rib(-15) > 0.3
      ; interpolate between levels 1 and 2
      y = (/Rib_mid(0), Rib_mid(1)/)
      x = (/-5., -15./)
      quad_slope = 0.
    else
      if (Rib_mid(2).gt.0.3) then
        ; Rib(-15) < 0.3, Rib(-25) > 0.3
        ; interpolate between levels 2 and 3
        y = (/Rib_mid(1), Rib_mid(2)/)
        x = (/-15., -25./)
        quad_slope = (Rib_mid(1) - Rib_mid(0))/-10.
      else
        ; Rib(-25) < 0.3, Rib(-35) > 0.3
        ; interpolate between levels 3 and 4
        y = (/Rib_mid(2), Rib_mid(3)/)
        x = (/-25., -35./)
        quad_slope = (Rib_mid(2) - Rib_mid(1))/-10.
      end if
    end if
    slope = (y(1) - y(0))/(x(1) - x(0))
    ; linear interpolation is easy
    ; (y(x)-y0) = slope*(x-x0)
    ; find x such that y(x) = 0.3:
    ; (0.3 - y0) = slope*(x-x0)
    ; x = x0+(0.3 - y0)/slope
    linear(i) = x(0)+(0.3-y(0))/slope
    if (extra_plot) then
      plot_Rib_mid = Rib_mid
      do ii=0,Ncells-2
        tmp_slope = (Rib_mid(ii+1)-Rib_mid(ii))/(midpoints(ii+1)-midpoints(ii))
;        print((/slope/))
        do jj=0, Nvert
          kk = ii*Nvert+jj
          lin_col(kk) = tmp_slope*(v_depth(kk)-midpoints(ii))+Rib_mid(ii)
          if (v_depth(kk).ge.hmix(i)) then
            exact(kk) = 0.
          else
            exact(kk) = 0.15*(hmix(i)-v_depth(kk))
          end if
;          print((/v_depth(kk)+": "+lin_col(kk)/))
        end do
      end do
    end if

    ; quadratic interpolation is tougher
    slope = quad_slope
    ; find a, b, c such that
    ; y = ax^2 + bx + c => y' = 2x+b
    ; y(x0)  = y0
    ; y(x1)  = y1
    ; y'(x0) = slope
    ; Then find x such that y(x) = 0.3 (x0 <= x <= x1)
    ; Matrix method
    ; [ x1^2 x1 1 ][ a ]   [    y1 ]
    ; [ x0^2 x0 1 ][ b ] = [    y0 ]
    ; [  2x0  1 0 ][ c ]   [ slope ]
    ; A little linear algebra yields
    ; [ a ]   [    1/(x1-x0)^2          -1/(x1-x0)^2    -1/(x1-x0)   ][    y1 ]
    ; [ b ] = [ -2x0/(x1-x0)^2         2x0/(x1-x0)^2 (x1+x0)/(x1-x0) ][    y0 ]
    ; [ c ]   [ (x0/(x1-x0))^2 -x1(2x0-x1)/(x1-x0)^2 -x1x0/(x1-x0)   ][ slope ]
    det = -((x(1)-x(0))^2)
    Minv = new((/3, 3/), float)
    col_vec = new((/3,1/), float)
    Minv(0,0) = -1.
    Minv(0,1) = 1.
    Minv(0,2) = x(1) - x(0)
    Minv(1,0) = 2.*x(0)
    Minv(1,1) = -2.*x(0)
    Minv(1,2) = x(0)^2 - x(1)^2
    Minv(2,0) = -(x(0)^2)
    Minv(2,1) = x(1)*(2.*x(0)-x(1))
    Minv(2,2) = x(1)*x(0)*(x(1)-x(0))
    col_vec = (/y(1), y(0), slope/)
    
    a = (/Minv(0,:)#col_vec/)/det
    b = (/Minv(1,:)#col_vec/)/det
    c = (/Minv(2,:)#col_vec/)/det
    delete(Minv)
    delete(col_vec)

    ; Error Checking
    print_lots = False
    if (abs(y(0) - (a*x(0)*x(0)+b*x(0)+c)).gt.tol) then
      print_lots = True
    end if
    if (abs(y(1) -(a*x(1)*x(1)+b*x(1)+c)).gt.tol) then
      print_lots = True
    end if
    if (abs(slope - (2*a*x(0)+b)).gt.tol) then
      print_lots = True
    end if
    if (print_lots) then
      print((/"MNL "+ y(0)+" "+(a*x(0)*x(0)+b*x(0)+c)/))
      print((/"MNL "+ y(1)+" "+(a*x(1)*x(1)+b*x(1)+c)/))
      print((/"MNL "+slope+" "+(2*a*x(0)+b)/))
    end if

    detM = b*b-4.*a*(c-0.3)
    if (detM.lt.-tol) then
      print((/"Error, no root to quadratic"/))
    end if
    if (abs(detM).le.tol) then
      quad(i) = -b/(2.*a)
    end if
    if (abs(detM).gt.tol) then
      root1 = (-b + sqrt(detM))/(2.*a)
      root2 = (-b - sqrt(detM))/(2.*a)
      ; Want x1 <= root <= x0
      num_roots = 0 ; "binary" counter for roots will be easy to extend to cubic
      if ((root1.gt.x(1)).and.(root1.lt.x(0))) then
        num_roots = num_roots+1
      end if
      if ((root2.gt.x(1)).and.(root2.lt.x(0))) then
        num_roots = num_roots+10
      end if

      if (num_roots.eq.0) then
        print((/"ERROR ("+detM+"): neither root ("+root1+", "+root2+"), is in ["+x(1)+","+x(0)+"]."/))
      else
        if (num_roots.eq.1) then
          quad(i) = root1
        else
          if (num_roots.eq.10) then
            quad(i) = root2
          else
            print((/"WARNING: two roots in ["+x(1)+","+x(0)+"] -- ("+root1+", "+root2+"), using first one."/))
          end if
        end if
      end if
    end if
    if (extra_plot) then
      do ii=0,Ncells-2
        tmp_slope = (Rib_mid(ii+1)-Rib_mid(ii))/(midpoints(ii+1)-midpoints(ii))
;        print((/slope/))
        do jj=0, Nvert
          kk = ii*Nvert+jj
          if ((Rib_mid(ii).lt.0.3).and.(Rib_mid(ii+1).ge.0.3)) then
            quad_col(kk) = v_depth(kk)*(a*v_depth(kk)+b)+c 
          else
            quad_col(kk) = tmp_slope*(v_depth(kk)-midpoints(ii))+Rib_mid(ii)
            if (v_depth(kk).ge.hmix(i)) then
              exact(kk) = 0.
            else
              exact(kk) = 0.15*(hmix(i)-v_depth(kk))
            end if
;            print((/v_depth(kk)+": "+lin_col(kk)/))
          end if
        end do
      end do
    end if

    ; Cubic splines are even more ridiculous
    slope2 = (y(1)-y(0))/(x(1)-x(0))
    ; find a, b, c, d such that
    ; y = ax^3 + bx^2 + cx + d => y' = 3ax^2 + 2bx + c
    ; y(x0)  = y0
    ; y(x1)  = y1
    ; y'(x0) = slope
    ; y'(x1) = slope2
    ; Then find x such that y(x) = 0.3 (x0 <= x <= x1)
    ; Matrix method
    ; [  x1^3 x1^2 x1 1 ][ a ]   [     y1 ]
    ; [  x0^3 x0^2 x0 1 ][ b ] = [     y0 ]
    ; [ 3x0^2  2x0  1 0 ][ c ]   [  slope ]
    ; [ 3x1^2  2x1  1 0 ][ d ]   [ slope2 ]
    ; A little linear algebra yields
    ; det = -(x1-x0)^3
; [a]   [             2            -2              x0-x1              x0-x1 ][     y1 ]
; [b] = [     -3(x1+x0)      3(x1+x0)    (x1-x0)(2x1+x0)    (x1-x0)(x1+2x0) ][     y0 ]
; [c]   [         6x1x0        -6x1x0 -x1(x1-x0)(x1+2x0) -x0(x1-x0)(2x1+x0) ][  slope ]
; [d]   [ -x0^2(3x1-x0) -x1^2(x1-3x0)      x0(x1-x0)x1^2      x0^2x1(x1-x0) ][ slope2 ]
    det = (x(0)-x(1))^3
    Minv = new((/4,4/), float)
    col_vec = new((/4,1/), float)
    Minv(0,0) = 2.
    Minv(0,1) = -2.
    Minv(0,2) = x(0)-x(1)
    Minv(0,3) = x(0)-x(1)

    Minv(1,0) = -3.*(x(1)+x(0))
    Minv(1,1) = 3.*(x(1)+x(0))
    Minv(1,2) = (x(1) - x(0))*(2.*x(1) + x(0))
    Minv(1,3) = (x(1) - x(0))*(x(1) + 2.*x(0))

    Minv(2,0) = 6.*x(1)*x(0)
    Minv(2,1) = -6.*x(1)*x(0)
    Minv(2,2) = -x(1)*(x(1) - x(0))*(x(1) + 2.*x(0))
    Minv(2,3) = -x(0)*(x(1) - x(0))*(2.*x(1) + x(0))

    Minv(3,0) = -(x(0)^2)*(3.*x(1)-x(0))
    Minv(3,1) = -(x(1)^2)*(x(1) - 3.*x(0))
    Minv(3,2) = x(1)^2*x(0)*(x(1)-x(0))
    Minv(3,3) = x(1)*x(0)^2*(x(1)-x(0))
    col_vec = (/y(1), y(0), slope, slope2/)

    a = (/Minv(0,:)#col_vec/)/det
    b = (/Minv(1,:)#col_vec/)/det
    c = (/Minv(2,:)#col_vec/)/det
    d = (/Minv(3,:)#col_vec/)/det
    delete(Minv)
    delete(col_vec)

    y1 = a*x(1)^3 + b*x(1)^2 + c*x(1) + d
    y0 = a*x(0)^3 + b*x(0)^2 + c*x(0) + d
    s0 = 3.*a*x(0)^2 + 2.*b*x(0) + c
    s1 = 3.*a*x(1)^2 + 2.*b*x(1) + c
    if (abs(y0-y(0)).gt.tol) then
      print((/"MNL y0: "+y(0)+" "+y0/))
    end if
    if (abs(y1-y(1)).gt.tol) then
      print((/"MNL y1: "+y(1)+" "+y1/))
    end if
    if (abs(s0-slope).gt.tol) then
        print((/"MNL slope: "+slope+" "+s0/))
    end if
    if (abs(s1-slope2).gt.tol) then
        print((/"MNL slope2: "+slope2+" "+s1/))
    end if

    ; Root find
    root = 0.5*(x(0) + x(1))
;    root = x(1)
    do j=1,100
      fval  = a*root^3 + b*root^2 + c*root + d-0.3
      if (abs(fval).lt.tol) then
        total_its = total_its + j
;        print((/"MNL MNL - Newton method converged in "+(j-1)+" iterations"/))
        break
      else
        sval = 3.*a*root^2 + 2.*b*root + c
        root = root - fval/sval
      end if
    end do
    cubic(i) = root
    if (extra_plot) then
      do ii=0,Ncells-2
        tmp_slope = (Rib_mid(ii+1)-Rib_mid(ii))/(midpoints(ii+1)-midpoints(ii))
;        print((/slope/))
        do jj=0, Nvert
          kk = ii*Nvert+jj
          if ((Rib_mid(ii).lt.0.3).and.(Rib_mid(ii+1).ge.0.3)) then
            cubic_col(kk) = v_depth(kk)*(v_depth(kk)*(a*v_depth(kk)+b)+c)+d 
          else
            cubic_col(kk) = tmp_slope*(v_depth(kk)-midpoints(ii))+Rib_mid(ii)
            if (v_depth(kk).ge.hmix(i)) then
              exact(kk) = 0.
            else
              exact(kk) = 0.15*(hmix(i)-v_depth(kk))
            end if
          end if
        end do
      end do
    end if


    print((/hmix(i)+" "+cont(i)+" "+linear(i)+" "+quad(i)+" "+cubic(i)/))
  end do
;  print((/"Total iterations: "+total_its/))

  wks = gsn_open_wks(out_type, "bldepth")
  res = True
  res@gsnMaximize = True
  ; Given gsnMaximize = True, HeightF and WidthF just set aspect ratio
  res@vpHeightF = 24
  res@vpWidthF = 20
  res@tiMainString = "Boundary Layer Depth"
  res@tiXAxisString = "HMIX (m)"
  res@tiYAxisString = "HBL (m)"
  res@xyLineThicknesses = (/0.5, 2., 1., 2./)
  res@xyDashPatterns    = (/  0,  0,  1,  1/)
  res@tmXBMode = "Explicit"
  res@tmYLMode = "Explicit"

  if (all_pos) then
    hmix = abs(hmix)
    cont = abs(cont)
    linear = abs(linear)
    quad = abs(quad)
    cubic = abs(cubic)
    res@trXMinF = 10.
    res@trXMaxF = 30.
    res@trYMinF = 8.
    res@trYMaxF = 32.
    res@tmXBValues = (/10, 15, 20, 25, 30/)
    res@tmXBLabels = (/10, 15, 20, 25, 30/)
    res@tmXBMinorValues = ispan(10, 30, 1)
    res@tmYLValues = (/10, 15, 20, 25, 30/)
    res@tmYLLabels = (/10, 15, 20, 25, 30/)
    res@tmYLMinorValues = ispan(8, 32, 1)
  else
    res@trXMinF = -30.
    res@trXMaxF = -10.
    res@trYMinF = -32.
    res@trYMaxF = -8.
    res@tmXBValues = (/-30, -25, -20, -15, -10/)
    res@tmXBLabels = (/-30, -25, -20, -15, -10/)
    res@tmXBMinorValues = ispan(-30, -10, 1)
    res@tmYLValues = (/-30, -25, -20, -15, -10/)
    res@tmYLLabels = (/-30, -25, -20, -15, -10/)
    res@tmYLMinorValues = ispan(-32, -8, 1)
  end if

  plot = gsn_csm_xy(wks, hmix, (/cont, linear, quad, cubic/), res)
  delete(wks)
  delete(res@tmYLValues)
  delete(res@tmYLLabels)
  delete(res@tmYLMinorValues)
  delete(res@vpWidthF)
  delete(res@vpHeightF)

  wks = gsn_open_wks(out_type, "bldepth_error")
  res@tiMainString = "Error in Boundary Layer Depth"
  res@tiYAxisString = "HBL - HMIX (m)"
  if (all_pos) then
    res@trYMinF = -0.6
    res@trYMaxF = 3.1
    res@tmYLValues = (/-0.5, 0., 0.5, 1., 1.5, 2., 2.5, 3./)
    res@tmYLLabels = (/-0.5, 0., 0.5, 1., 1.5, 2., 2.5, 3./)
    res@tmYLMinorValues = ispan(-6, 31, 1)/10.
  else
    res@trYMinF = -3.1
    res@trYMaxF = 0.6
    res@tmYLValues = (/-3., -2.5, -2., -1.5, -1., -0.5,  0., 0.5/)
    res@tmYLLabels = (/-3., -2.5, -2., -1.5, -1., -0.5,  0., 0.5/)
    res@tmYLMinorValues = ispan(-31, 6, 1)/10.
  end if
  
  plot = gsn_csm_xy(wks, hmix, (/cont-cont, cont-linear, cont-quad, cont-cubic/), res)
  delete(wks)
  delete(res)

  if (make_last_plot) then
    wks = gsn_open_wks(out_type, "single_col"+abs(plot_at))
    res = True
    res@gsnDraw = False
    res@gsnFrame = False

    res@tiMainString = "Bulk Richardson number"
    fmt="%6.2f"
    res@gsnCenterString = "                " + \
                          "(hmix = "+sprintf(fmt, plot_at)+" m, hbl = "+     \
                          sprintf(fmt, plot_at-2)+" m):C:" +                 \
                          "linear = "+sprintf(fmt, -abs(linear(icol))) +     \
                          " m, quad = "+sprintf(fmt, -abs(quad(icol))) +     \
                          " m , cubic = " + sprintf(fmt, -abs(cubic(icol)))+ \
                          " m"
    res@tiXAxisString = "Ri_b"
    res@tiYAxisString = "depth (m)"
    res@trXMinF = 0.
    res@trXMaxF = plot_Rib_mid(Ncells-1)
    res@xyLineThicknesses = (/0.5, 0.5, 2., 1., 2./)
    res@xyDashPatterns    = (/  2,   2,  0,  1,  1/)
    res@gsnMaximize = True
  
    plot  = gsn_csm_xy(wks, (/Ri_crit, exact, lin_col, quad_col, cubic_col/), v_depth, res)

    res2 = res
    delete(res2@xyLineThicknesses)
    delete(res2@xyDashPatterns)
    res2@xyMarkLineMode = "Markers"
    res2@xyMonoMarker = True
    res2@xyMarker = 16
    plot2 = gsn_csm_xy(wks, plot_Rib_mid, midpoints, res2)
    overlay(plot, plot2)
    draw(plot)
    frame(wks)
  end if

  print((/"Best:"/))
  imin = min(ind(abs(linear-cont).eq.min(abs(linear-cont))))
  print((/hmix(imin)+": "+abs(linear(imin)-cont(imin))/))
  imin = min(ind(abs(quad-cont).eq.min(abs(quad-cont))))
  print((/hmix(imin)+": "+abs(quad(imin)-cont(imin))/))
  imin = min(ind(abs(cubic-cont).eq.min(abs(cubic-cont))))
  print((/hmix(imin)+": "+abs(cubic(imin)-cont(imin))/))

  print((/"Worst:"/))
  imin = min(ind(abs(linear-cont).eq.max(abs(linear-cont))))
  print((/hmix(imin)+": "+abs(linear(imin)-cont(imin))/))
  imin = min(ind(abs(quad-cont).eq.max(abs(quad-cont))))
  print((/hmix(imin)+": "+abs(quad(imin)-cont(imin))/))
  imin = min(ind(abs(cubic-cont).eq.max(abs(cubic-cont))))
  print((/hmix(imin)+": "+abs(cubic(imin)-cont(imin))/))

end
