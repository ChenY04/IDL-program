function mapinfo,name
  envi_open_file, name, r_fid=fid, /no_realize
  envi_file_query, fid, ns = ns, nl = nl, nb = nb, Data_Type = Data_Type, dims=dims
  map_info = envi_get_map_info(fid=fid)
  return,map_info
end

function get_landsat_date,dir,y_p,d_p
  dir_l = dir+strtrim(y_p,1)+'\*'
  file_l = file_search(dir_l)
  n_la = n_elements(file_l) & L = STRARR(n_la/2)
  for di = 0, n_la/2-1 do begin
    L[di] = file_l[di*2]
  endfor
  n_l = strmid(L,49,51) & n1 = fix(n_l)
  d = abs(n1-d_p) & min_d = min(d)
  return,min_d
end

;========================all prediction main
pro all_data_prediction ,y_p,d_p,fsdaf_weight,ci,predictionname
  ;t0=systime(1)

  year = [2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016]

  w_r = fltarr(400,400,16)
  w_all = fltarr(400,400,16)
  temp_wr = fltarr(400,400,16)
  in_year_w = fltarr(400,400,16)
  image_e = fltarr(400,400,16)
  image_a = fltarr(400,400,16)
  weight_fix = fltarr(400,400,16)
  k_l=fltarr(400,400,16) & lc_w = fltarr(400,400,16)
  in_year_w_temp = fltarr(400,400,16)
  ns = 400 & nl = 400
  image_f = fltarr(400,400,16)
  cloud_weight = fltarr(400,400)
  L_c = fltarr(400,400) & new_coe = fltarr(400,400,16)


  map_info = mapinfo('G:\final_data\evergreen_ndvi\modis\2001\mod09a1_2001_1')
  dir = 'G:\final_data\evergreen_ndvi\landsat\'
  dir_fsdaf = 'G:\final_data\evergreen_ndvi\FSDAF_other_other\'
  ;===============================================================
  for yeari = 0,15 do begin
    
    w_r_name = fsdaf_weight+'fsdaf_weight_r_'+strtrim(year[yeari],1)+'_'+strtrim(d_p,1)+'_'+strtrim(ci,1)
    envi_open_file, w_r_name, r_fid=fid, /no_realize
    envi_file_query, fid,dims=dims
    w_r[*,*,yeari] = abs(envi_get_data(fid=fid,dims=dims,pos=0))

    kl_name = 'G:\final_data\evergreen_ndvi\change_landsat\0\mod09a1_change_'+strtrim(year[yeari],1)+'_'+strtrim(d_p,1)+'_'+strtrim(ci,1)
    envi_open_file, kl_name, r_fid=fid, /no_realize
    envi_file_query, fid,dims=dims
    k_l[*,*,yeari] = envi_get_data(fid=fid,dims=dims,pos=0)

  endfor

  ;===============================in_year
  print, 'get the in_year weight'

  for i = 0,15 do begin
    temp_wr[*,*,i]=w_r[*,*,i]/total(w_r,3)
  endfor
  a = mean(temp_wr)
  for i = 0,ns-1 do begin
    for j = 0,nl-1 do begin
      for b = 0,16-1 do begin
        in_year_w_temp[i,j,b] = a^temp_wr[i,j,b]
      endfor
    endfor
  endfor
  for i = 0,15 do begin
    in_year_w[*,*,i]=in_year_w_temp[*,*,i]/total(in_year_w_temp,3)
  endfor

  ;===================================

  envi_write_envi_file, in_year_w, out_name='G:\final_data\evergreen_ndvi\final_weight\enhance\in_year_weight' $
  +strtrim(y_p,1)+strtrim(d_p,1), data_type=4,$
    ns=400, nl=400, nb=16, map_info = map_info, /no_open

  ;===================
  print, 'get the landcover_change weight'

  coe = fltarr(ns,nl,16)
  coe = landc_weight(y_p,d_P)
  
  for i = 0,15 do begin
    coe_temp = fltarr(400,400)
    coe_temp = coe[*,*,i]
    L_cname = 'G:\final_data\evergreen_ndvi\classification\0\mod09a1_class_'+strtrim(year[i],1)+'_'+strtrim(d_p,1)
    envi_open_file, L_cname, r_fid=fid, /no_realize
    envi_file_query, fid,dims=dims
    L_c[*,*] = envi_get_data(fid=fid,dims=dims,pos=0)
    new_coe[*,*,i] = relation_landsat(L_c,coe_temp)
  endfor

  for i = 0,15 do begin
    lc_w[*,*,i]=new_coe[*,*,i]/total(new_coe,3)
  endfor
  ;==========================


  envi_write_envi_file, lc_w, out_name='G:\final_data\evergreen_ndvi\final_weight\enhance\landcoverchange_weight' $
  +strtrim(y_p,1)+strtrim(d_p,1), data_type=4,$
    ns=400, nl=400, nb=16, map_info = map_info, /no_open

  print,'change weight finished'
  ;===================
  min_d = get_landsat_date(dir,y_p,d_p)
  temp_w_all =  in_year_w*lc_w
  if min_d le 8 then begin
    weight_fix[*,*,y_p-2001] = 1.0
    temp_w_all = weight_fix
  endif
  ;======================

  for i = 0,15 do begin
    w_all[*,*,i]=temp_w_all[*,*,i]/total(temp_w_all,3)
  endfor

  envi_write_envi_file, w_all, out_name='G:\final_data\evergreen_ndvi\final_weight\enhance\prediction_weight_' $
  +strtrim(y_p,1)+strtrim(d_p,1), data_type=4,$
    ns=400, nl=400, nb=16, map_info = map_info, /no_open


  for i = 0,15 do begin
    name_f = dir_fsdaf+'mod09a1_'+strtrim(year[i],1)+'_'+strtrim(d_p,1)+'_FSDAF_'+strtrim(ci,1)
    envi_open_file, name_f, r_fid=fid, /no_realize
    envi_file_query, ns=ns,nl=nl,fid,dims=dims
    image_f[*,*,i] = envi_get_data(fid=fid,dims=dims,pos=0)
  endfor


  for i = 0,15 do begin
    temp_f = image_f[*,*,i]
    temp_kl = fltarr(400,400)
    if min_d gt 8 then begin
    if abs(mean(k_l[*,*,i])) lt 0.20 then begin
      temp_kl = k_l[*,*,i]
    endif
    endif
    cloud_weight = cloud_qa(d_p,i)
    image_e[*,*,i] = (temp_f+temp_kl)*cloud_weight
    image_a[*,*,i] = image_e[*,*,i]*w_all[*,*,i]
  endfor

  ;========================================================
  final_prediction = total(image_a,3)

  envi_write_envi_file, final_prediction, out_name=predictionname, data_type=4,$
    ns=400, nl=400, nb=1, map_info = map_info, /no_open

  print, 'all the data prediction finished'
  ;print, 'all data process time used:', floor((systime(1)-t0)/3600), 'h',floor(((systime(1)-t0) mod 3600)/60),'m',(systime(1)-t0) mod 60,'s'

END
