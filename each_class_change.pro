Pro GetData,ImgData = ImgData,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type,$
  FileName = FileName,Map_info = map_Info, Fid = Fid
  Filter = ['all file;*.*']
  Envi_Open_File,FileName,R_Fid = Fid
  Envi_File_Query,Fid,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type
  map_info = envi_get_map_info(fid=Fid)
  dims = [-1,0,ns - 1 ,0,nl - 1]
  case Data_Type Of
    1:ImgData = BytArr(ns,nl,nb)
    2:ImgData = IntArr(ns,nl,nb)
    3:ImgData = LonArr(ns,nl,nb)
    4:ImgData = FltArr(ns,nl,nb)
    5:ImgData = DblArr(ns,nl,nb)
    6:ImgData = COMPLEXARR(ns,nl,nb)
    9:ImgData = DCOMPLEXARR(ns,nl,nb)
    12:ImgData = UINTARR(ns,nl,nb)
    13:ImgData = ULONARR(ns,nl,nb)
    14:ImgData = LON64ARR(ns,nl,nb)
    15:ImgData = ULON64ARR(ns,nl,nb)
  EndCase
  For i = 0,nb-1 Do Begin
    Dt = Envi_Get_Data(Fid = Fid,dims = dims,pos=i)
    ImgData[*,*,i] = Dt[*,*]
  EndFor
End

;function of multiple linear regression without interception
Function P_OLS, ind_v,dep_v,k,low,high    ;format: ind_v:k*n,dep_v:1*n

  common V_PUB, matrix
  common V_PUB2, y

  nb=n_elements(dep_v)
  x=reform(dep_v)
  varend=fltarr(k)
  yend=fltarr(k,nb)
  for ii=0,k-1 do begin
    varend[ii]=sqrt(total((ind_v[ii,*]-mean(ind_v[ii,*]))^2))
    yend[ii,*]=ind_v[ii,*]
  endfor
  var=sqrt(total((dep_v-mean(dep_v))^2))

  y=dep_v
  matrix=yend
  y=transpose(double(y))

  glow=fltarr(k)+low
  ghigh=fltarr(k)+high
  gintial=fltarr(k)+1.0/float(k)
  gbnd   = [glow, ghigh]
  Lbnd =[0,100]
  nobj   = 1
  g      = gintial
  Lcomp = 'HMBL11'
  nobj=0
  CONSTRAINED_MIN, g, gbnd, Lbnd, nobj, Lcomp, inform
  L =total((matrix ## g- y)^2)
  return,g
END

FUNCTION HMBL11, g
  common V_PUB
  common V_PUB2
  L=total((matrix ## g-y)^2)
  RETURN, L
END


;-------------------------------------------------------------------
;                       main program
;-------------------------------------------------------------------

function each_class_change,L1_name,M1_name,M2_name,year,d_p

  ;parameters
  ;----------------------------------------------------------------------
  min_class=2.0
  max_class=3.0
  num_pure=100
  DN_min=-1.0
  DN_max=1.0
  scale_factor=16
  block_size=30
  background=-2
  background_band=1
  temp_file='E:\temp'
  ;------------------------------------------------------------------------
  ;define the variables
  ;-----------------------------------
  k_l=fltarr(400,400)
  L_c = fltarr(400,400)

  ;-------------------------------------------------------------------------------------

    ;open the fine image of the first pair
    FileName1 = L1_name
    GetData,ImgData=fine1_whole,FileName = FileName1,Fid=fid0
    ;envi_open_file,FileName1,r_fid=fid0
    envi_file_query,fid0,ns=ns,nl=nl,nb=nb,dims=dims


    map_info = envi_get_map_info(fid = fid0)
    patch_long=block_size*scale_factor
    orig_ns=ns
    orig_nl=nl
    n_ns=ceil(float(orig_ns)/patch_long)
    n_nl=ceil(float(orig_nl)/patch_long)

    ind_patch1=intarr(4,n_ns*n_nl)           ;divide the whole scene into 1000*1000 block
    ind_patch=intarr(4,n_ns*n_nl)
    location=intarr(4,n_ns*n_nl)

    for i_ns=0,n_ns-1,1 do begin
      for i_nl=0,n_nl-1,1 do begin
        ind_patch1[0,n_ns*i_nl+i_ns]=i_ns*patch_long
        ind_patch[0,n_ns*i_nl+i_ns]=max([0,ind_patch1[0,n_ns*i_nl+i_ns]-scale_factor])
        location[0,n_ns*i_nl+i_ns]=ind_patch1[0,n_ns*i_nl+i_ns]-ind_patch[0,n_ns*i_nl+i_ns]

        ind_patch1[1,n_ns*i_nl+i_ns]=min([ns-1,(i_ns+1)*patch_long-1])
        ind_patch[1,n_ns*i_nl+i_ns]=min([ns-1,ind_patch1[1,n_ns*i_nl+i_ns]+scale_factor])
        location[1,n_ns*i_nl+i_ns]=ind_patch1[1,n_ns*i_nl+i_ns]-ind_patch1[0,n_ns*i_nl+i_ns]+location[0,n_ns*i_nl+i_ns]

        ind_patch1[2,n_ns*i_nl+i_ns]=i_nl*patch_long
        ind_patch[2,n_ns*i_nl+i_ns]=max([0,ind_patch1[2,n_ns*i_nl+i_ns]-scale_factor])
        location[2,n_ns*i_nl+i_ns]=ind_patch1[2,n_ns*i_nl+i_ns]-ind_patch[2,n_ns*i_nl+i_ns]

        ind_patch1[3,n_ns*i_nl+i_ns]=min([nl-1,(i_nl+1)*patch_long-1])
        ind_patch[3,n_ns*i_nl+i_ns]=min([nl-1,ind_patch1[3,n_ns*i_nl+i_ns]+scale_factor])
        location[3,n_ns*i_nl+i_ns]=ind_patch1[3,n_ns*i_nl+i_ns]-ind_patch1[2,n_ns*i_nl+i_ns]+location[2,n_ns*i_nl+i_ns]
      endfor
    endfor

    tempoutname=temp_file+'\temp_F1'

    pos=indgen(nb)
    for isub=0,n_ns*n_nl-1,1 do begin
      dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
      envi_doit, 'resize_doit', fid=fid0, pos=pos, dims=dims, interp=0, rfact=[1,1], $
        out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
      envi_file_mng, id=r_fid1, /remove
    endfor




    ;open the coarse image of the first pair
    ;-----------------------------------------------------------
    FileName2 = M1_name
    envi_open_file,FileName2,r_fid=fid
    tempoutname=temp_file+'\temp_C1'
    pos=indgen(nb)
    for isub=0,n_ns*n_nl-1,1 do begin
      dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
      envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
        out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
      envi_file_mng, id=r_fid1, /remove
    endfor
    envi_file_mng, id=fid, /remove


    ;open the coarse image of the prediction time
    ;-----------------------------------------------------------
    FileName5 = M2_name
    envi_open_file,FileName5,r_fid=fid
    tempoutname=temp_file+'\temp_C0'
    pos=indgen(nb)
    for isub=0,n_ns*n_nl-1,1 do begin
      dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
      envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
        out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
      envi_file_mng, id=r_fid1, /remove
    endfor
    envi_file_mng, id=fid, /remove



    ;replace background pixels by mean of non-background to avoid its effects on classification
    background_whole=bytarr(ns,nl)
    ind_back=where(fine1_whole[*,*,background_band-1] eq background, num_back)
    if (num_back gt 0) then begin
      background_whole[ind_back]=1
      for iband=0, nb-1, 1 do begin
        temp=fine1_whole[*,*,iband]
        temp[ind_back]=mean(temp[where(background_whole eq 0)])
        fine1_whole[*,*,iband]=temp
      endfor
    endif
    tempoutname11=temp_file+'\fine1_nobackground'
    Envi_Write_Envi_File,fine1_whole,Out_Name = tempoutname11
    ind_back=0
    temp=0
    fine1_whole=0
    background_whole=0; clear variable

    envi_open_file,tempoutname11,r_fid=fid00


    ;step1: get spectral classes from fine resolution image at t1 by isodata
    ;parameter of isodata
    CHANGE_THRESH = .05
    NUM_CLASSES = max_class
    ITERATIONS = 20
    ISO_MERGE_DIST = 0.05*DN_max
    ISO_MERGE_PAIRS = 2
    ISO_MIN_PIXELS = 200
    ISO_SPLIT_SMULT = 1
    ISO_SPLIT_STD = 0.05*DN_max
    MIN_CLASSES = min_class
    out_bname = 'IsoData'
    out_name=temp_file+'\class_ISODATA'
    ENVI_File_Query, fid00, DIMS=dims, NB=nb1
    ENVI_DOIT, 'class_doit', fid=fid00, pos=indgen(nb1), dims=dims, $
      out_bname=out_bname, out_name=out_name, method=4, $
      r_fid=r_fid, $
      NUM_CLASSES = NUM_CLASSES, $
      ITERATIONS = ITERATIONS, $
      CHANGE_THRESH = CHANGE_THRESH, $
      ISO_MERGE_DIST = ISO_MERGE_DIST, $
      ISO_MERGE_PAIRS = ISO_MERGE_PAIRS, $
      ISO_MIN_PIXELS = ISO_MIN_PIXELS, $
      ISO_SPLIT_SMULT = ISO_SPLIT_SMULT, $
      ISO_SPLIT_STD = ISO_SPLIT_STD, $
      MIN_CLASSES = MIN_CLASSES


    envi_open_file,out_name,r_fid=fid
    envi_file_query, fid, dims=dims
    ;ENVI_DOIT, 'CLASS_MAJORITY_DOIT', CLASS_PTR=array, DIMS=dims, FID=fid, /IN_MEMORY, KERNEL_SIZE=array, METHOD=0, OUT_BNAME='majority',  POS=[0], R_FID=fid
    
    ;=========================================================
    L_c[*,*] = envi_get_data(fid=fid,dims=dims,pos=0)
    envi_write_envi_file, L_c, out_name='G:\final_data\evergreen_ndvi\classification\0\mod09a1_class_'+strtrim(year,1)+'_'+strtrim(d_p,1), data_type=4,$
      ns=400, nl=400, nb=1, /no_open
    

    tempoutname=temp_file+'\class'
    pos=indgen(1)
    for isub=0,n_ns*n_nl-1,1 do begin
      dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
      envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
        out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
      envi_file_mng, id=r_fid1, /remove
    endfor

    envi_file_mng, id=fid, /remove,/delete
    envi_file_mng, id=fid00, /remove,/delete
    envi_file_mng, id=fid0, /remove
    ;

    ;------------------------------------------------------------------
    ;process  each block
    ;-------------------------------------------------------------------


    ;print,'there are total',n_ns*n_nl,' blocks'

    for isub=0,n_ns*n_nl-1,1 do begin

      ;open each block image *****

      FileName = temp_file+'\temp_F1'
      GetData,ImgData=fine1,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type,FileName = FileName+strtrim(isub+1,1),Fid = Fid11
      fine1=float(fine1)

      FileName = temp_file+'\temp_C1'
      GetData,ImgData=coarse1,FileName = FileName+strtrim(isub+1,1),Fid = Fid12
      coarse1=FLOAT(coarse1)

      FileName = temp_file+'\temp_C0'
      GetData,ImgData=coarse2,FileName = FileName+strtrim(isub+1,1),Fid = Fid13
      coarse2=FLOAT(coarse2)

      FileName = temp_file+'\class'
      GetData,ImgData=L1_class0,FileName = FileName+strtrim(isub+1,1),Fid = Fid14

      num_class=max(L1_class0)
      ;recode the classification map if the subset does not have all classes
      i_new_c=0
      L1_class=intarr(ns,nl)
      for iclass=0, num_class-1,1 do begin
        ind_ic=where(L1_class0 eq iclass+1 and fine1[*,*, background_band-1] ne background, num_ic)
        if (num_ic gt 0) then begin
          L1_class[ind_ic]=i_new_c+1
          i_new_c=i_new_c+1
        endif
      endfor

      num_class=max(L1_class)


      if (num_class gt 0) then begin   ;do not process if the whole subset is background  ****


        ;correct extreme noise in fine1 becase extreme values will affect the allowed data range
        for ib=0,nb-1, 1 do begin
          sortIndex = Sort(fine1[*,*,ib])
          sortIndices = (Findgen(float(ns)*nl+1))/(float(ns)*nl)
          Percentiles=[0.0001, 0.9999]
          dataIndices = Value_Locate(sortIndices, Percentiles)
          data_1_4= (fine1[*,*,ib])[sortIndex[dataIndices]]
          ;correct too small values
          ind_small=where(fine1[*,*,ib] le data_1_4[0] or fine1[*,*,ib] lt DN_min)
          temp=fine1[*,*,ib]
          temp[ind_small]=min((fine1[*,*,ib])[where(fine1[*,*,ib] gt data_1_4[0] and fine1[*,*,ib] ge DN_min)])
          fine1[*,*,ib]=temp
          ;correct too large values
          ind_large=where(fine1[*,*,ib] ge data_1_4[1] or fine1[*,*,ib] gt DN_max)
          temp=fine1[*,*,ib]
          temp[ind_large]=max((fine1[*,*,ib])[where(fine1[*,*,ib] lt data_1_4[1] and fine1[*,*,ib] le DN_max)])
          fine1[*,*,ib]=temp
        endfor


        ;get index image between coarse and fine resolutions
        ii=0
        ns_c=floor(ns/scale_factor)
        nl_c=floor(nl/scale_factor)
        index_f=intarr(ns,nl)
        index_c=intarr(ns_c,nl_c)
        for i=0, ns_c-1, 1 do begin
          for j=0,nl_c-1,1 do begin
            index_f[i*scale_factor:(i+1)*scale_factor-1, j*scale_factor:(j+1)*scale_factor-1]=ii
            index_c[i,j]=ii
            ii=ii+1.0
          endfor
        endfor

        ;col and row index
        row_ind=intarr(ns,nl)
        col_ind=intarr(ns,nl)
        for i=0,ns-1,1 do begin
          col_ind[i,*]=i
        endfor
        for i=0,nl-1,1 do begin
          row_ind[*,i]=i
        endfor

        ;resample coarse image to coarse resolution
        fine_c1=fltarr(ns_c,nl_c,nb)
        coarse_c1=fltarr(ns_c,nl_c,nb)
        coarse_c2=fltarr(ns_c,nl_c,nb)
        row_c=fltarr(ns_c,nl_c)
        col_c=fltarr(ns_c,nl_c)
        for ic=0,ns_c-1, 1 do begin
          for jc=0,nl_c-1, 1 do begin
            ind_c=where(index_f eq index_c[ic,jc])
            row_c[ic,jc]= mean(row_ind[ind_c])
            col_c[ic,jc]= mean(col_ind[ind_c])
            for ib=0,nb-1,1 do begin
              fine_c1[ic,jc,ib]=mean((fine1[*,*,ib])[ind_c])
              coarse_c1[ic,jc,ib]=mean((coarse1[*,*,ib])[ind_c])
              coarse_c2[ic,jc,ib]=mean((coarse2[*,*,ib])[ind_c])
            endfor
          endfor
        endfor


        ;print,'number of class:',num_class

        ;step 2: get fracture of each class within each coarse pixel at t1
        Fraction1=fltarr(ns_c,nl_c,num_class)
        for ic=0,ns_c-1, 1 do begin
          for jc=0,nl_c-1, 1 do begin
            ind_c=where(index_f eq index_c[ic,jc], num_c)
            L1_class_c=L1_class[ind_c]
            for iclass=0, num_class-1,1 do begin
              ind_ic=where(L1_class_c eq iclass+1, num_ic)
              Fraction1[ic,jc,iclass]=float(num_ic)/float(num_c)
            endfor
            if (total(Fraction1[ic,jc,*]) le 0.999) then begin   ;avoild pixels have background fine pixels
              Fraction1[ic,jc,*]=0
            endif
          endfor
        endfor


        ;step 3: METHOD2:estimate average spectral change of each class using pixels without land cover change
        c_rate=fltarr(num_class,nb)

        ;allowed change value for each band
        min_allow=fltarr(nb)
        max_allow=fltarr(nb)
        for ib=0,nb-1,1 do begin
          min_allow[ib]=min(coarse_c2[*,*,ib]-coarse_c1[*,*,ib])-stddev(coarse_c2[*,*,ib]-coarse_c1[*,*,ib])
          max_allow[ib]=max(coarse_c2[*,*,ib]-coarse_c1[*,*,ib])+stddev(coarse_c2[*,*,ib]-coarse_c1[*,*,ib])
        endfor

        for ib=0,nb-1, 1 do begin
          X_matrix=fltarr(num_class,num_pure*num_class)
          y_matrix=fltarr(1,num_pure*num_class)
          ii=0
          for ic=1, num_class, 1 do begin
            order=sort(Fraction1[*,*,ic-1])
            order=reverse(order)
            ind_f=where(Fraction1[*,*,ic-1] gt 0.01, num_f) ;make sure all selected modis pixel contain class i
            num_pure1=min([num_f,num_pure])
            change_c=(coarse_c2[*,*,ib])[order[0:(num_pure1-1)]]-(coarse_c1[*,*,ib])[order[0:(num_pure1-1)]]

            ;only use 0.1-0.9 samples to exclude the land cover change pixels
            sortIndex = Sort(change_c)
            sortIndices = (Findgen(num_pure1+1))/float(num_pure1)
            Percentiles=[0.25, 0.75]
            ;Percentiles=[0.1, 0.9]
            dataIndices = Value_Locate(sortIndices, Percentiles)
            data_1_4= change_c[sortIndex[dataIndices]]
            ind_nonchange=where(change_c ge data_1_4[0] and change_c le data_1_4[1],num_nonc)
            y_matrix[0,ii:ii+num_nonc-1]=change_c[ind_nonchange]
            for icc=0, num_class-1,1 do begin
              f_c=(Fraction1[*,*,icc])[order[0:(num_pure1-1)]]
              X_matrix[icc,ii:ii+num_nonc-1]=f_c[ind_nonchange]
            endfor
            ii=ii+num_nonc
          endfor
          X_matrix=X_matrix[*,0:ii-1]
          y_matrix=y_matrix[0,0:ii-1]
          result=P_OLS(X_matrix,y_matrix,num_class,min_allow[ib],max_allow[ib])  ;format: ind_v:k*n,dep_v:1*n
          c_rate[*,ib]=result[*]
        endfor
      endif
    endfor

    temp_kl = fltarr(ns,nl)
    for bi=0,ns-1 do begin
      for bj=0,nl-1 do begin
        temp_kl[bi,bj]=c_rate[L_c[bi,bj]-1]
      endfor
    endfor
    k_l=temp_kl
    
    return,k_l

end
    