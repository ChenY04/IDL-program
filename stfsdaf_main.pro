pro STFSDAF_main
  t1=systime(1)

  ;y_p = 2007
  year = [2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016]
  dir_m = 'G:\final_data\evergreen_ndvi\modis\'
  dir = 'G:\final_data\evergreen_ndvi\landsat\'
  fsdaf_weight = 'G:\final_data\evergreen_ndvi\fsdaf_weight\'
  ;predictionname = 'G:\final_data\evergreen_ndvi\result\mod09a1_'+strtrim(y_P,1)+strtrim(d_p,1)+'_STFSDAF_v32'

  ;==========================================
  nm1 = indgen(46)*8+1
  n= n_elements(nm1)

;  for bi= 0,13 do begin
;    print,'processing the date:'
;    d_p = nm1[bi]
;    print,d_p
;
;    for yeari = 0,15 do begin
;      print,yeari
;      print,'doing fsdaf:'
;      dir_l = dir+strtrim(year[yeari],1)+'\*'
;      file_l = file_search(dir_l)
;      fsdaf_weight_name_r = fsdaf_weight+'fsdaf_weight_r_'+strtrim(year[yeari],1)+'_'+strtrim(d_p,1)+'_'+strtrim(bi,1)
;      n_la = n_elements(file_l)
;      L = STRARR(n_la/2)
;      for i = 0, n_la/2-1 do begin
;        L[i] = file_l[i*2]
;      endfor
;      n_l = strmid(L,49,51)
;      n1 = fix(n_l)
;      d = abs(n1-d_p)
;      min_d = min(d,index)
;      l_date = n1[index]
;      File_Name1 = L[index]
;      nm = indgen(46)*8+1
;      d2 = abs(nm-l_date)
;      min_d2 = min(d2,index2)
;      m_date = nm[index2]
;      print,L[index]
;      print,m_date
;      File_Name2 = dir_m+strtrim(year[yeari],1)+'\mod09a1_'+strtrim(year[yeari],1)+'_'+strtrim(m_date,1)
;      File_Name5 = dir_m+strtrim(year[yeari],1)+'\mod09a1_'+strtrim(year[yeari],1)+'_'+strtrim(d_p,1)
;      batch_fsdaf,File_Name1,File_Name2,File_Name5,year[yeari],fsdaf_weight_name_r,d_p,bi
;    endfor
;    endfor


    for yeari =2,2 do begin
    y_p = year[yeari]
    for ci = 0,45 do begin
    print,nm1[ci]
    d_p = nm1[ci]
    print,'doing class change'
    landc_chang_landsat,y_p,d_p,ci
    print,'doing all prediction:'
    predictionname = 'G:\final_data\evergreen_ndvi\result\mod09a1_'+strtrim(y_p,1)+strtrim(d_p,1)+'_STFSDAF_v1'
    all_data_prediction,y_p,d_p,fsdaf_weight,ci,predictionname
    endfor
    endfor

  print,'===================='
  print, 'all the program finished time used:', floor((systime(1)-t1)/3600), 'h',floor(((systime(1)-t1) mod 3600)/60),'m',(systime(1)-t1) mod 60,'s'


end