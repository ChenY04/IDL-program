function cloud_qa,d_p,yeari
 year = [2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016]
 dir_cloud = 'G:\final_data\evergreen_ndvi\cloud_qa\'
 dir = 'G:\final_data\evergreen_ndvi\landsat\'
 
 dir_l = dir+strtrim(year[yeari],1)+'\*'
 file_l = file_search(dir_l)
 n_la = n_elements(file_l)
 L = STRARR(n_la/2)
 for i = 0, n_la/2-1 do begin
   L[i] = file_l[i*2]
 endfor
 n_l = strmid(L,49,51)
 l_name = strmid(L,42,51)
 n1 = fix(n_l)
 d = abs(n1-d_p)
 min_d = min(d,index)
 cloud_weight = fltarr(400,400)
 envi_open_file, dir_cloud+l_name[index]+'_qa', r_fid=fid, /no_realize
 envi_file_query, ns=ns,nl=nl,fid,dims=dims
 cloud_weight[*,*] = envi_get_data(fid=fid,dims=dims,pos=0)
 a=fltarr(400,400)+1.0
 cloud_weight = a-cloud_weight
 
 return,cloud_weight
end
