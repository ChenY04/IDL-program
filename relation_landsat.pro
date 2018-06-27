function relation_landsat,L_c,lc_temp

  scale_factor = 16
  ns = 400 & nl = 400
  
  ns_c=floor(ns/scale_factor)
  nl_c=floor(nl/scale_factor)
  index_f=fltarr(ns,nl)
  index_r=fltarr(ns,nl)
  temp_class=fltarr(ns_c,nl_c)
  num_class = max(L_c)
  for i=0, ns_c-1, 1 do begin
    for j=0,nl_c-1,1 do begin
      temp_class = L_c[i*scale_factor:((i+1)*scale_factor-1), j*scale_factor:((j+1)*scale_factor-1)]
      for k =1.0,num_class do begin
      rate_pixel = float(n_elements(where(temp_class eq k)))/float(n_elements(temp_class))
      if rate_pixel gt 0.9 then begin
        index_f[i*scale_factor:(i+1)*scale_factor-1, j*scale_factor:(j+1)*scale_factor-1] = k
        index_r[i*scale_factor:(i+1)*scale_factor-1, j*scale_factor:(j+1)*scale_factor-1] = lc_temp[i*scale_factor+scale_factor/2,j*scale_factor+scale_factor/2]
      endif
      endfor
    endfor
  endfor
  for i = 1.0,num_class do begin
  L_c[where(L_c eq i)]=mean(index_r[where(index_f eq i)])
  endfor
  return,L_c
end