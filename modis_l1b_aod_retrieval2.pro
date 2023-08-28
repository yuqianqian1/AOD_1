function hdf4_data_get,file_name,sds_name
  sd_id=hdf_sd_start(file_name,/read)
  sds_index=hdf_sd_nametoindex(sd_id,sds_name)
  sds_id=hdf_sd_select(sd_id,sds_index)
  hdf_sd_getdata,sds_id,data
  hdf_sd_endaccess,sds_id
  hdf_sd_end,sd_id
  return,data
end

function hdf4_caldata_get,file_name,sds_name,scale_name,offset_name
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
  data_cal=fltarr(data_size[1],data_size[2],data_size[3])
  for layer_i=0,data_size[3]-1 do data_cal[*,*,layer_i]=scale_data[layer_i]*(data[*,*,layer_i]-offset_data[layer_i])
  data=!null
  return,data_cal
end

function aod_lookup,lut_sz,lut_vz,lut_ra,lut_aod,col_num,line_num,lut_data,sz_angle,vz_angle,ra_angle,band_ref,surf_band,cloud_result,fr_ref
  sz_n=n_elements(lut_sz)
  vz_n=n_elements(lut_vz) 
  ra_n=n_elements(lut_ra)
  aod_n=n_elements(lut_aod)
  Combination_i=vz_n*ra_n
  Combination_j=sz_n*vz_n*ra_n
  band_aod=fltarr(col_num,line_num)
  for line_i=0,line_num-1 do begin
    for col_i=0,col_num-1 do begin
      if (cloud_result[col_i,line_i] eq 1) or (fr_ref[col_i,line_i] lt 0.01) or (fr_ref[col_i,line_i] gt 0.25) then continue
      aod_lut_subset=fltarr(4,aod_n)
      lut_band_tri=fltarr(7,8)

      sz_pos_stat=sort(abs(sz_angle[col_i,line_i]-lut_sz))
      sz_pos0=lut_sz[sz_pos_stat[0]]
      sz_pos1=lut_sz[sz_pos_stat[1]]

      vz_pos_stat=sort(abs(vz_angle[col_i,line_i]-lut_vz))
      vz_pos0=lut_vz[vz_pos_stat[0]]
      vz_pos1=lut_vz[vz_pos_stat[1]]

      ra_pos_stat=sort(abs(ra_angle[col_i,line_i]-lut_ra))
      ra_pos0=lut_ra[ra_pos_stat[0]]
      ra_pos1=lut_ra[ra_pos_stat[1]]

      for aod_i=0,aod_n-1 do begin
        comb_pos=0
        for sz_pos=0,1 do begin
          for vz_pos=0,1 do begin
            for ra_pos=0,1 do begin
              lut_band_tri[*,comb_pos]=lut_data[*,(sz_pos_stat[sz_pos]*combination_i)+(vz_pos_stat[vz_pos]*ra_n)+ra_pos_stat[ra_pos]+(ulong64(aod_i)*combination_j)]
              comb_pos=comb_pos+1
            endfor
          endfor
        endfor

        ;对观测几何三维插值
        for prm_i=0,2 do begin
          xd=float(sz_angle[col_i,line_i]-sz_pos0)/float(sz_pos1-sz_pos0)
          yd=float(vz_angle[col_i,line_i]-vz_pos0)/float(vz_pos1-vz_pos0)
          zd=float(ra_angle[col_i,line_i]-ra_pos0)/float(ra_pos1-ra_pos0)
          aod_lut_subset[prm_i,aod_i]=lut_band_tri[prm_i,0]*(1.0-xd)*(1.0-yd)*(1.0-zd)+$
            lut_band_tri[prm_i,1]*(1.0-xd)*(1.0-yd)*(zd)+$
            lut_band_tri[prm_i,2]*(1.0-xd)*(yd)*(1.0-zd)+$
            lut_band_tri[prm_i,3]*(1.0-xd)*(yd)*(zd)+$
            lut_band_tri[prm_i,4]*(xd)*(1.0-yd)*(1.0-zd)+$
            lut_band_tri[prm_i,5]*(xd)*(1.0-yd)*(zd)+$
            lut_band_tri[prm_i,6]*(xd)*(yd)*(1.0-zd)+$
            lut_band_tri[prm_i,7]*(xd)*(yd)*(zd)
        endfor
        aod_lut_subset[3,aod_i]=lut_aod[aod_i]
      endfor

      toa_sim=aod_lut_subset[0,*]+(aod_lut_subset[1,*]*surf_band[col_i,line_i])/(1.0-surf_band[col_i,line_i]*aod_lut_subset[2,*])
      delta_toa=toa_sim-band_ref[col_i,line_i]

      gt0_pos=where(delta_toa gt 0)
      lt0_pos=where(delta_toa lt 0)
      if (gt0_pos[0] eq -1) or (lt0_pos[0] eq -1) then continue
      max_neg_delta=max(delta_toa[lt0_pos])
      min_posi_delta=min(delta_toa[gt0_pos])
      max_neg_pos=where(delta_toa eq max_neg_delta)
      min_posi_pos=where(delta_toa eq min_posi_delta)
      max_neg_aod=aod_lut_subset[3,max_neg_pos]
      min_posi_aod=aod_lut_subset[3,min_posi_pos]
      k=(max_neg_aod-min_posi_aod)/(max_neg_delta-min_posi_delta)
      aod_interp=k*(0.0-min_posi_delta)+min_posi_aod
      band_aod[col_i,line_i]=aod_interp
    endfor
  endfor
  return,band_aod
end

pro data_glt,out_res,data_col,data_line,target_lon_data,target_lat_data,target_data,data_box_geo_out,geo_info 
  congrid_scale=8
  target_lon_data_interp=congrid(target_lon_data,data_col*congrid_scale,data_line*congrid_scale,/interp)
  target_lat_data_interp=congrid(target_lat_data,data_col*congrid_scale,data_line*congrid_scale,/interp)

  lon_min=min(target_lon_data)
  lon_max=max(target_lon_data)
  lat_min=min(target_lat_data)
  lat_max=max(target_lat_data)

  data_box_geo_col=ceil((lon_max-lon_min)/out_res)
  data_box_geo_line=ceil((lat_max-lat_min)/out_res)
  data_box_lon=fltarr(data_box_geo_col,data_box_geo_line)
  data_box_lon[*,*]=-9999.0
  data_box_lat=fltarr(data_box_geo_col,data_box_geo_line)
  data_box_lat[*,*]=-9999.0
  data_box_geo=fltarr(data_box_geo_col,data_box_geo_line)

  data_box_lon_col_pos=floor((target_lon_data_interp-lon_min)/out_res)
  data_box_lon_line_pos=floor((lat_max-target_lat_data_interp)/out_res)
  data_box_lon[data_box_lon_col_pos,data_box_lon_line_pos]=target_lon_data_interp
  
  data_box_lat_col_pos=floor((target_lon_data_interp-lon_min)/out_res)
  data_box_lat_line_pos=floor((lat_max-target_lat_data_interp)/out_res)
  data_box_lat[data_box_lon_col_pos,data_box_lon_line_pos]=target_lat_data_interp
  
  data_box_geo_col_pos=floor((target_lon_data-lon_min)/out_res)
  data_box_geo_line_pos=floor((lat_max-target_lat_data)/out_res)
  data_box_geo[data_box_geo_col_pos,data_box_geo_line_pos]=(target_data gt 0.0)*target_data+(target_data le 0.0)*(-9999.0)

  data_box_geo_out=fltarr(data_box_geo_col,data_box_geo_line)
  window_size=9
  jump_size=(window_size-1)/2
  for data_box_geo_line_i=jump_size,data_box_geo_line-jump_size-1 do begin
    for data_box_geo_col_i=jump_size,data_box_geo_col-jump_size-1 do begin
      if data_box_geo[data_box_geo_col_i,data_box_geo_line_i] eq 0.0 then begin
        distance=sqrt((data_box_lon[data_box_geo_col_i,data_box_geo_line_i]-data_box_lon[(data_box_geo_col_i-jump_size):(data_box_geo_col_i+jump_size),(data_box_geo_line_i-jump_size):(data_box_geo_line_i+jump_size)])^2+$
          (data_box_lat[data_box_geo_col_i,data_box_geo_line_i]-data_box_lat[(data_box_geo_col_i-jump_size):(data_box_geo_col_i+jump_size),(data_box_geo_line_i-jump_size):(data_box_geo_line_i+jump_size)])^2)
        distance_sort_pos=sort(distance)
        data_box_geo_window=data_box_geo[(data_box_geo_col_i-jump_size):(data_box_geo_col_i+jump_size),(data_box_geo_line_i-jump_size):(data_box_geo_line_i+jump_size)]
        data_box_geo_sort=data_box_geo_window[distance_sort_pos]
        fill_pos=where(data_box_geo_sort ne 0.0)
        fill_value=data_box_geo_sort[fill_pos[0]]
        data_box_geo_out[data_box_geo_col_i,data_box_geo_line_i]=fill_value
      endif else begin
        data_box_geo_out[data_box_geo_col_i,data_box_geo_line_i]=data_box_geo[data_box_geo_col_i,data_box_geo_line_i]
      endelse
    endfor
  endfor

  data_box_geo_out=abs((data_box_geo_out gt 0.0)*data_box_geo_out)*(data_box_lat ne -9999.0)

  geo_info={$
    MODELPIXELSCALETAG:[out_res,out_res,0.0],$
    MODELTIEPOINTTAG:[0.0,0.0,0.0,lon_min,lat_max,0.0],$
    GTMODELTYPEGEOKEY:2,$
    GTRASTERTYPEGEOKEY:1,$
    GEOGRAPHICTYPEGEOKEY:4326,$
    GEOGCITATIONGEOKEY:'GCS_WGS_1984',$
    GEOGANGULARUNITSGEOKEY:9102,$
    GEOGSEMIMAJORAXISGEOKEY:6378137.0,$
    GEOGINVFLATTENINGGEOKEY:298.25722}
end

pro modis_l1b_aod_retrieval
  start_time=systime(1)
  modis_file='P:/coarse_data/chapter_5/MYD021KM.A2019116.0630.061.2019117022052.hdf'
  cloud_file='P:/coarse_data/chapter_5/MYD35_L2.A2019116.0630.061.2019117022432.hdf'
  result_name='P:/coarse_data/chapter_5/MYD021KM.A2019116.0630.061.2019117022052_aod_geo.tiff'
  output_directory='P:/coarse_data/chapter_5/'
  red_lut_file='P:/coarse_data/chapter_5/modis_lut_ref_t_red_omar4_mls.txt'
  blue_lut_file='P:/coarse_data/chapter_5/modis_lut_ref_t_blue_omar4_mls.txt'
  qkm_ref=hdf4_caldata_get(modis_file,'EV_250_Aggr1km_RefSB','reflectance_scales','reflectance_offsets')
  hkm_ref=hdf4_caldata_get(modis_file,'EV_500_Aggr1km_RefSB','reflectance_scales','reflectance_offsets')
  km_ems=hdf4_caldata_get(modis_file,'EV_1KM_Emissive','radiance_scales','radiance_offsets')
  ref_size=size(qkm_ref)
  
  modis_lon_data=congrid(hdf4_data_get(modis_file,'Longitude'),ref_size[1],ref_size[2],/interp)
  modis_lat_data=congrid(hdf4_data_get(modis_file,'Latitude'),ref_size[1],ref_size[2],/interp)
  
  sz_angle=0.01*hdf4_data_get(modis_file,'SolarZenith')
  sa_angle=0.01*hdf4_data_get(modis_file,'SolarAzimuth')
  vz_angle=0.01*hdf4_data_get(modis_file,'SensorZenith')
  va_angle=0.01*hdf4_data_get(modis_file,'SensorAzimuth')
  ra_angle=abs(sa_angle-va_angle)
  ra_angle=(ra_angle le 180.0)*ra_angle+(ra_angle gt 180.0)*(360.0-ra_angle)

  sz_angle=congrid(sz_angle,ref_size[1],ref_size[2],/interp)
  vz_angle=congrid(vz_angle,ref_size[1],ref_size[2],/interp)
  ra_angle=congrid(ra_angle,ref_size[1],ref_size[2],/interp)

  scatter_angle_cos=-(cos(vz_angle*!DTOR)*cos(sz_angle*!DTOR)+sin(vz_angle*!DTOR)*sin(sz_angle*!DTOR)*cos(ra_angle*!DTOR))
  scatter_angle=acos(scatter_angle_cos)/!DTOR
  
  red_ref=qkm_ref[*,*,0]/cos(sz_angle*!DTOR)
  blue_ref=hkm_ref[*,*,0]/cos(sz_angle*!DTOR)
  nr_ref=hkm_ref[*,*,2]/cos(sz_angle*!DTOR)
  fr_ref=hkm_ref[*,*,4]/cos(sz_angle*!DTOR)
  qkm_ref=!null
  hkm_ref=!null
  
  cloud_data=hdf4_data_get(cloud_file,'Cloud_Mask')
  cloud_0=cloud_data[*,*,0]
  cloud_0=(cloud_0 ge 0)*cloud_0+(cloud_0 lt 0)*(128+abs(cloud_0))
  cloud_0_size=size(cloud_0)
  
  cloud_binary=bytarr(cloud_0_size[1],cloud_0_size[2],8)
  for cloud_i=0,7 do begin
    cloud_binary[*,*,cloud_i]=cloud_0 mod 2
    cloud_0=cloud_0/2
  endfor
  
  cloud_result=(cloud_binary[*,*,0] eq 1) and (cloud_binary[*,*,1] eq 0) and (cloud_binary[*,*,2] eq 0)
  
  ndvi_swir=(nr_ref-fr_ref)/(nr_ref+fr_ref)
  slope_ndvi=(ndvi_swir lt 0.25)*0.48+$
    (ndvi_swir gt 0.75)*0.58+$
    ((ndvi_swir ge 0.25) and (ndvi_swir le 0.75))*(0.48+0.2*(ndvi_swir-0.25))
  
  slope=slope_ndvi+0.002*scatter_angle-0.27
  yint=-0.00025*scatter_angle+0.033
  
  surf_red=slope*fr_ref+yint
  surf_blue=0.49*surf_red+0.005
  
  openr,1,red_lut_file
  red_lut_data=fltarr(7,file_lines(red_lut_file))
  readf,1,red_lut_data
  free_lun,1
  
  openr,1,blue_lut_file
  blue_lut_data=fltarr(7,file_lines(blue_lut_file))
  readf,1,blue_lut_data
  free_lun,1
  
  sz=[0,6,12,18,24,30,36,42,48,54,60,66,72,78]
  vz=[0,6,12,18,24,30,36,42,48,54,60,66,72,78]
  ra=[0,12,24,36,48,60,72,84,96,108,120,132,144,156,168,180]
  aod=[0.01,0.25,0.5,1.0,1.25,1.5,2.0,3.0,5.0]

  red_aod=aod_lookup(sz,vz,ra,aod,ref_size[1],ref_size[2],red_lut_data,sz_angle,vz_angle,ra_angle,red_ref,surf_red,cloud_result,fr_ref)
  blue_aod=aod_lookup(sz,vz,ra,aod,ref_size[1],ref_size[2],blue_lut_data,sz_angle,vz_angle,ra_angle,blue_ref,surf_blue,cloud_result,fr_ref)
  modis_target_data=(red_aod+blue_aod)/2.0
  resolution=0.01
  data_glt,resolution,ref_size[1],ref_size[2],modis_lon_data,modis_lat_data,modis_target_data,data_box_geo_out,geo_info 

  print,result_name
  write_tiff,result_name,data_box_geo_out,/float,geotiff=geo_info
  end_time=systime(1)
  print,'Time consuming: '+strtrim(string(end_time-start_time),2)+' s.'
end