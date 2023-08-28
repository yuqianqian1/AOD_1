function hdf4_caldata_get,file_name,sds_name,scale_name,offset_name
;重新封装后的定标函数
  sd_id=hdf_sd_start(file_name,/read)
  sds_index=hdf_sd_nametoindex(sd_id,sds_name)
  sds_id=hdf_sd_select(sd_id,sds_index)
  hdf_sd_getdata,sds_id,data
  
  att_index=hdf_sd_attrfind(sds_id,scale_name)
  hdf_sd_attrinfo,sds_id,att_index,data=scale_data
  att_index=hdf_sd_attrfind(sds_id,offset_name)
  hdf_sd_attrinfo,sds_id,att_index,data=offset_data
  
  hdf_sd_endaccess,sds_id
  hdf_sd_end,sd_id
  
  data_size=size(data)
  data_ref=fltarr(data_size[1],data_size[2],data_size[3])
  for layer_i=0,data_size[3]-1 do begin
    data_ref[*,*,layer_i]=scale_data[layer_i]*(data[*,*,layer_i]-offset_data[layer_i])
  endfor
  data=!null
  return,data_ref
end

function hdf4_data_get,file_name,sds_name
  sd_id=hdf_sd_start(file_name,/read)
  sds_index=hdf_sd_nametoindex(sd_id,sds_name)
  sds_id=hdf_sd_select(sd_id,sds_index)
  hdf_sd_getdata,sds_id,data
  hdf_sd_endaccess,sds_id
  hdf_sd_end,sd_id
  return,data
end

;主函数
pro modis_l1b_aod_retrieval
  compile_opt idl2
  envi,/restore_base_save_files
  envi_batch_init
  
  modis_file='F:/论文/气溶胶/数据下载/data/2月28号/MYD021KM.A2022059.0550.061.2022059165345.hdf'
  result_name='F:/论文/气溶胶/数据下载/data/2月28号/MYD021KM.A2022059.0550.061.2022059165345_AOD_GEO.hdf'
  ;云腌摸产品数据的获取
  cloud_file='F:/论文/气溶胶/数据下载/data/2月28号/MYD35_L2.A2022059.0550.061.2022059165621.hdf'
  output_directory='F:/论文/气溶胶/数据下载/data/2月28号/result/'
  red_lut_file='F:/论文/气溶胶/数据下载/data/2月28号/modis_lut_red.txt'
  blue_lut_file='F:/论文/气溶胶/数据下载/data/2月28号/modis_lut_blue.txt'
  km_ref=hdf4_caldata_get(modis_file,'EV_1KM_RefSB','reflectance_scales','reflectance_offsets')
  qkm_ref=hdf4_caldata_get(modis_file,'EV_250_Aggr1km_RefSB','reflectance_scales','reflectance_offsets')
  hkm_ref=hdf4_caldata_get(modis_file,'EV_500_Aggr1km_RefSB','reflectance_scales','reflectance_offsets')
  ;注意这里的代码和以前调用的数据不一样了
  km_ems=hdf4_caldata_get(modis_file,'EV_1KM_Emissive','radiance_scales','radiance_offsets')
  ref_size=size(qkm_ref)
  
  
  ;对各种角做插值
  sz_angle=0.01*hdf4_data_get(modis_file,'SolarZenith')
  sa_angle=0.01*hdf4_data_get(modis_file,'SolarAzimuth')
  vz_angle=0.01*hdf4_data_get(modis_file,'SensorZenith')
  va_angle=0.01*hdf4_data_get(modis_file,'SensorAzimuth')
  
  ra_angle=abs(sa_angle-va_angle)
  ra_angle=(ra_angle le 180.0)*ra_angle+(ra_angle gt 180.0)*(360.0-ra_angle)
  
  sz_angle=congrid(sz_angle,ref_size[1],ref_size[2],/interp)
  vz_angle=congrid(vz_angle,ref_size[1],ref_size[2],/interp)
  ;sa_angle=congrid(hdf4_data_get(modis_file,'SolarAzimuth')*0.01,ref_size[1],ref_size[2],/interp)
  ra_angle=congrid(ra_angle,ref_size[1],ref_size[2],/interp)
  
  
  sca_angle_cos=-(cos(vz_angle*!DTOR)*cos(sz_angle*!DTOR)+sin(vz_angle*!DTOR)*sin(sz_angle*!DTOR)*cos(ra_angle*!DTOR))
  sca_angle=acos(sca_angle_cos)/!DTOR
  
  
  
  ;气溶胶反演所需数据
  red_ref=qkm_ref[*,*,0]/cos(sz_angle*!DTOR)
  blue_ref=hkm_ref[*,*,0]/cos(sz_angle*!DTOR)
  nr_ref=hkm_ref[*,*,2]/cos(sz_angle*!DTOR)
  fr_ref=hkm_ref[*,*,4]/cos(sz_angle*!DTOR)
  ;获取数据集大小,任意一个波段的行列数是一样的
  
  
  qkm_ref=!null
  hkm_ref=!null
  modis_lon_data=hdf4_data_get(modis_file,'Longitude')
  modis_lat_data=hdf4_data_get(modis_file,'Latitude')
  
  ;云腌摸数据的提取
  cloud_data=hdf4_data_get(cloud_file,'Cloud_Mask')
  ;转成二进制
  cloud_0=cloud_data[*,*,0]
  cloud_0=(cloud_0 ge 0)*cloud_0+(cloud_0 lt 0)*(128+abs(cloud_0))
  cloud_0_size=size(cloud_0)
  
  cloud_binary=bytarr(cloud_0_size[1],cloud_0_size[2],8)
  for cloud_i=0,7 do begin
    cloud_binary[*,*,cloud_i]=cloud_0 mod 2
    cloud_0=cloud_0/2
  endfor
  
  ;满足清晰无云的条件的数据
  
 clear_result=(cloud_binary[*,*,0]eq 1)and(cloud_binary[*,*,1]eq 0)and(cloud_binary[*,*,2]eq 0)
 ;波段运算
  ndvi_swir=(nr_ref-fr_ref)/(nr_ref+fr_ref)
  slope=(ndvi_swir lt 0.25)*0.48+$
    (ndvi_swir ge 0.75)*0.58+$
    ((ndvi_swir ge 0.25)and(ndvi_swir le 0.75))*(0.48+0.2*(ndvi_swir-0.25))
  

  
  slope=slope+0.002*sca_angle-0.27
  yint=-0.00025*sca_angle+0.033
  ;红波段的地表反射率
  surf_red=fr_ref*slope+yint
  surf_blue=surf_red*0.49+0.005
  
  ;查找近似行
  openr,1,red_lut_file
  red_lut_data=fltarr(7,file_lines(red_lut_file))
  readf,1,red_lut_data
  free_lun,1
  
  openr,1,blue_lut_file
  blue_lut_data=fltarr(7,file_lines(blue_lut_file))
  readf,1,blue_lut_data
  free_lun,1
  
  sz=[0,10,20,30,40,50,60,70,80]
  sz_n=n_elements(sz)
  vz=[0,10,20,30,40,50,60,70,80]
  vz_n=n_elements(vz)
  ra=[0,12,24,36,48,60,72,84,96,108,120,132,144,156,168,180]
  ra_n=n_elements(ra)
  aod=[0.01,0.25,0.5,1.0,1.25,1.5,2.0,3.0,5.0]
  aod_n=n_elements(aod)

  red_aod=fltarr(ref_size[1],ref_size[2])
  blue_aod=fltarr(ref_size[1],ref_size[2])
  ;print,red_lut_data[*,500]
  ;print,blue_lut_data[*,500]
  ;red
  for line_i=0,ref_size[2]-1 do begin
    for col_i=0,ref_size[1]-1 do begin
      if (clear_result[col_i,line_i] eq 1)or(nr_ref[col_i,line_i] gt 0.25) or (nr_ref[col_i,line_i] lt 0.01) then continue
      ;插值处理
      aod_lut_subset=fltarr(7,aod_n)

      sz_temp=sz_angle[col_i,line_i]
      vz_temp=vz_angle[col_i,line_i]
      ra_temp=ra_angle[col_i,line_i]
                  
      ra_dif=abs(ra_temp-ra)     
      vz_dif=abs(vz_temp-vz) 
      sz_dif=abs(sz_temp-sz)
      
         
      ra_pos=where(ra_dif eq min(ra_dif))
      vz_pos=where(vz_dif eq min(vz_dif))
      sz_pos=where(sz_dif eq min(sz_dif))
      ;help,ra_pos ,它本身是long 
     
      for aod_i=0,aod_n-1 do begin
         line_pos=ra_pos[0]+vz_pos[0]*ra_n+sz_pos[0]*ra_n*vz_n+aod_i*ra_n*sz_n*vz_n
         aod_lut_subset[*,aod_i]=red_lut_data[*,line_pos]
      endfor
      ;print,sz_temp,vz_temp,ra_temp
      ;print,aod_lut_subset
      ;print,'-*****************-'
      
      toa_sim=aod_lut_subset[0,*]+(aod_lut_subset[1,*]*surf_red[col_i,line_i])/(1.0-aod_lut_subset[2,*]*surf_red[col_i,line_i])
      delta_toa=toa_sim-red_ref[col_i,line_i]
      
      gt0_pos=where(delta_toa gt 0)
      lt0_pos=where(delta_toa lt 0)
      if(gt0_pos[0] eq -1)or(lt0_pos[0] eq -1)then continue
      max_neg_delta=max(delta_toa[lt0_pos])
      min_posi_delta=min(delta_toa[gt0_pos])
      
      max_neg_pos=where(delta_toa eq max_neg_delta)
      min_posi_pos=where(delta_toa eq min_posi_delta)
      
      max_neg_aod=aod_lut_subset[6,max_neg_pos]
      min_posi_aod=aod_lut_subset[6,min_posi_pos]
      
      fit_k=(max_neg_aod-min_posi_aod)/(max_neg_delta-min_posi_delta)
      aod_interp=fit_k*(0.0-max_neg_delta)+max_neg_aod
      red_aod[col_i,line_i]=aod_interp
      ;2
      ;print,max_neg_pos,min_posi_pos
      ;print,'********'
      ;1
      ;print,max_neg_delta,min_posi_delta
      ;print,'********'
      
    endfor
  endfor
  red_aod=(red_aod gt 0)*red_aod
  ;write_tiff,output_directory+'aod_red2.tiff',red_aod,/float
  ;got重投影
  ;实时获取
  target_data_size=size(red_aod)
  modis_lon_data=congrid(modis_lon_data,target_data_size[1],target_data_size[2],/interp)
  modis_lat_data=congrid(modis_lat_data,target_data_size[1],target_data_size[2],/interp)

  out_lon=output_directory+'lon_out.tiff'
  out_lat=output_directory+'lat_out.tiff'
  out_target=output_directory+'target_red_aod.tiff'
  write_tiff,out_lon,modis_lon_data,/float
  write_tiff,out_lat,modis_lat_data,/float
  write_tiff,out_target,red_aod,/float

  envi_open_file,out_lon,r_fid=x_fid;打开经度数据，获取经度文件id,lon_out
  envi_open_file,out_lat,r_fid=y_fid;打开纬度数据，获取经度文件id,lat_out
  envi_open_file,out_target,r_fid=target_fid;打开目标数据，获取目标文件id

  out_name_glt=output_directory+file_basename(modis_file,'.hdf')+'_glt.img'
  out_name_glt_hdr=output_directory+file_basename(modis_file,'.hdf')+'_glt.hdr'
  
  i_proj=envi_proj_create(/geographic)
  o_proj=envi_proj_create(/geographic)
  envi_glt_doit,$
    i_proj=i_proj,x_fid=x_fid,y_fid=y_fid,x_pos=0,y_pos=0,$;指定创建GLT所需输入数据信息
    o_proj=o_proj,pixel_size=pixel_size,rotation=0.0,out_name=out_name_glt,r_fid=glt_fid;指定输出GLT文件信息

  out_name_geo=output_directory+file_basename(modis_file,'.hdf')+'_georef.img'
  out_name_geo_hdr=output_directory+file_basename(modis_file,'.hdf')+'_georef.hdr'
  envi_georef_from_glt_doit,$
    glt_fid=glt_fid,$;指定重投影所需GLT文件信息
    fid=target_fid,pos=0,$;指定待投影数据id
    out_name=out_name_geo,r_fid=geo_fid;指定输出重投影文件信息

  envi_file_query,geo_fid,dims=data_dims
  target_data=envi_get_data(fid=geo_fid,pos=0,dims=data_dims)

  map_info=envi_get_map_info(fid=geo_fid)
  geo_loc=map_info.(1)
  px_size=map_info.(2)

  geo_info={$
    MODELPIXELSCALETAG:[px_size[0],px_size[1],0.0],$
    MODELTIEPOINTTAG:[0.0,0.0,0.0,geo_loc[2],geo_loc[3],0.0],$
    GTMODELTYPEGEOKEY:2,$
    GTRASTERTYPEGEOKEY:1,$
    GEOGRAPHICTYPEGEOKEY:4326,$
    GEOGCITATIONGEOKEY:'GCS_WGS_1984',$
    GEOGANGULARUNITSGEOKEY:9102,$
    GEOGSEMIMAJORAXISGEOKEY:6378137.0,$
    GEOGINVFLATTENINGGEOKEY:298.25722}

  write_tiff,result_name,target_data,/float,geotiff=geo_info

  envi_file_mng,id=x_fid,/remove
  envi_file_mng,id=y_fid,/remove
  envi_file_mng,id=target_fid,/remove
  envi_file_mng,id=glt_fid,/remove
  envi_file_mng,id=geo_fid,/remove
  file_delete,[out_lon,out_lat,out_target,out_name_glt,out_name_glt_hdr,out_name_geo,out_name_geo_hdr]
  envi_batch_exit,/no_confirm
    
end