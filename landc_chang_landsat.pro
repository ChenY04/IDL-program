pro landc_chang_landsat,y_p,d_p,bi

  year = [2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016]
  ;------------------------------------------------------------------------

  ;file_path
  ;---------------------------------------------------------------------------------
  dir_m = 'G:\final_data\evergreen_ndvi\modis\'
  dir = 'G:\final_data\evergreen_ndvi\landsat\'
  k_l = fltarr(400,400,16)
  ;---------------------------------------------------------------------------------
    
  for yeari = 0,15 do begin
    print,year[yeari]
    
    dir_l = dir+strtrim(year[yeari],1)+'\*'
    file_l = file_search(dir_l)
    n_la = n_elements(file_l) & L = STRARR(n_la/2)

    for di = 0, n_la/2-1 do begin
      L[di] = file_l[di*2]
    endfor

    n_l = strmid(L,49,51) & n1 = fix(n_l)
    d = abs(n1-d_p) & min_d = min(d,index)
    L1_name = L[index];the fine image of the first pair
    M1_name = dir_m+strtrim(year[yeari],1)+'\mod09a1_'+strtrim(year[yeari],1)+'_'+strtrim(d_p,1);the coarse image of the first pair
    M2_name = dir_m+strtrim(y_p,1)+'\mod09a1_'+strtrim(y_p,1)+'_'+strtrim(d_p,1);the coarse image of the prediction time
    k_l[*,*,yeari] = each_class_change(L1_name,M1_name,M2_name,year[yeari],d_p);M2-M1
  endfor

  
  for i = 0,15 do begin
    envi_write_envi_file, k_l[*,*,i], out_name='G:\final_data\evergreen_ndvi\change_landsat\0\mod09a1_change_' $
    +strtrim(2001+i,1)+'_'+strtrim(d_p,1)+'_'+strtrim(bi,1), data_type=4,$
    ns=400,nl=400, nb=1, /no_open
  endfor
  
  end

