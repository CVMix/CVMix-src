; Basic settings that are the same for all reg_test plot features

GetFuncCode = ":"

undef("CVMix_res")
function CVMix_res(xmin,xmax,ymin,ymax)
begin
  return_val          = True
  return_val@gsnDraw  = False
  return_val@gsnFrame = False

  return_val@tiMainFuncCode  = GetFuncCode
  return_val@tiXAxisFuncCode = return_val@tiMainFuncCode
  return_val@tiYAxisFuncCode = return_val@tiMainFuncCode

  return_val@trXMinF = xmin
  return_val@trXMaxF = xmax
  return_val@trYMinF = ymin
  return_val@trYMaxF = ymax

  return(return_val)
end

undef("CVMix_txres")
function CVMix_txres(fontheight)
begin
  return_val = True
  return_val@txFontHeightF = fontheight
  return_val@txFuncCode    = GetFuncCode

  return(return_val)
end
