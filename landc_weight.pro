function landc_weight,y_p,d_P

year = [2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016]
dir_p = 'G:\final_data\evergreen_ndvi\time_series\mod09a1_'+strtrim(y_p,1)


envi_open_file, dir_p, r_fid=fid, /no_realize
envi_file_query, fid, ns = ns, nl = nl, nb = nb, Data_Type = Data_Type, dims=dims
img = fltarr(ns,nl,nb)
coe = fltarr(ns,nl,16)
for i = 0,nb-1 do begin
img[*,*,i] =envi_get_data(fid = fid,dims = dims,pos=i)
endfor

for yeari = 0,15 do begin
dir_o = 'G:\final_data\evergreen_ndvi\time_series\mod09a1_'+strtrim(year[yeari],1)
envi_open_file, dir_o, r_fid=fid, /no_realize
envi_file_query, fid, ns = ns, nl = nl, nb = nb, Data_Type = Data_Type, dims=dims
img_o = fltarr(ns,nl,nb)
for i = 0,nb-1 do begin
  img_o[*,*,i] =envi_get_data(fid = fid,dims = dims,pos=i)
endfor

point_p = (d_p-1)/8
b1 = point_p-3 & b2 = point_p+3
if b1 lt 0 then begin
  b1 = 0
endif
if b2 gt 45 then begin
  b2 = 45
endif

for i = 0,ns-1 do begin
  ;print,'sample:',i
for j = 0,nl-1 do begin

vector_p = reform(img[i,j,*])
vector_o = reform(img_o[i,j,*])
vector_p1 = vector_p[b1:b2] & vector_o1 = vector_o[b1:b2]
index1 = vector_p1 ge 0 and vector_p1 le 1
index2 = vector_o1 ge 0 and vector_o1 le 1
index = index1*index2
vector_p2 = vector_p1*index & vector_o2 = vector_o1*index
vector_p2 = vector_p2(where(vector_p2 ne 0)) & vector_o2 = vector_o2(where(vector_o2 ne 0))
;====================================

result = CORRELATE(vector_p2, vector_o2)
coe_temp = result[0]
coe[i,j,yeari]=coe_temp

endfor
endfor
endfor

for i = 0,ns-1 do begin
  for j = 0,nl-1 do begin
    num = coe[i,j,*]
    index = num ge 0
    num=num*index
    pos = where(~finite(num))
    pos_index = min(pos)
    if pos_index eq -1 then begin
      num = num
    endif else begin
    num[pos]=1 
    endelse
    coe[i,j,*] = num
  endfor  
endfor


envi_write_envi_file, coe, out_name='G:\final_data\evergreen_ndvi\final_weight\enhance\coe_'+strtrim(y_p,1)+strtrim(d_p,1), data_type=4,$
  ns=400, nl=400, nb=16, /no_open
  
return,coe

end