## This produces climatology data frames and matched BGC-argo data


library(data.table)
library(RCurl)
library(ggplot2)
library(ncdf4)
library(cwhmisc)
library(stringr)


dirData = "/perm_storage/home/kbaldry/Lennart/Data"
BGCArgo = fread(file.path(dirData, "SOCCOM_MLR_60S.csv"), header = T, sep = ",")
BGCArgo$lon[BGCArgo$lon>180] = BGCArgo$lon[BGCArgo$lon>180] - 360

url = "http://data.remss.com/ccmp/v02.0"



BGCArgo$dd = as.numeric(substr(BGCArgo$datetime_UTC,9,10))
date_list = paste(BGCArgo$yr,sprintf("%02.0f",BGCArgo$mm), sprintf("%02.0f",BGCArgo$dd),sep = "")
u_date_list = unique(date_list)
dir = "/perm_storage/home/kbaldry/Lennart/Data/Wind_data"
u_date_list = u_date_list[-which(as.numeric(u_date_list)>20190401)]
u_date_list = u_date_list[is.finite(as.numeric(u_date_list))]
BGCArgo$ws = NA
BGCArgo$sst = NA
BGCArgo$sea_ice = NA

for(x in u_date_list)
{
  # download files
  # wind speed
  file_location = file.path(url, paste("Y", substr(x,1,4), sep = ""),paste("M", substr(x,5,6), sep = ""))
  download.file(file.path(file_location, paste("CCMP_Wind_Analysis_",x,"_V02.0_L3.0_RSS.nc", sep = "")), file.path(dir, paste("CCMP_Wind_Analysis_",x,"_V02.0_L3.0_RSS.nc", sep = "")))
  
  
  
  # read files
  fl1 = nc_open(file.path(dir, paste("CCMP_Wind_Analysis_",x,"_V02.0_L3.0_RSS.nc", sep = "")))
  
  # create raster with mean daily windspeeds
  # should we based on time?
  avWS = raster(extent(range(ncvar_get(fl1,"longitude")) + c(-0.125,+0.125),range(ncvar_get(fl1,"latitude"))+ c(-0.125,+0.125)), resolution = 0.25) 
  avWS = setValues(avWS,t(apply(sqrt((ncvar_get(fl1,"uwnd")[,,1:4])^2+(ncvar_get(fl1,"vwnd")[,,1:4])^2),c(1,2),mean)))
  
  # create a raster with SST
  SST = readsst(as.Date(x, format = "%Y%m%d") ,time.resolution = "daily")
  SI = readice_daily(as.Date(x, format = "%Y%m%d"),setNA = F)
  AP = readwind(as.Date(x, format = "%Y%m%d"), time.resolution = "daily")
  
  # find index of mathing dates
  idx = which(date_list == x)
  # Bilinear interpolation
  # wind speed - 0-360 lon
  lon_ws = BGCArgo$lon[idx]
  lon_ws[lon_ws<=0] = lon_ws[lon_ws<=0] + 360
  BGCArgo$ws[idx] = raster::extract(avWS,cbind(lon_ws, BGCArgo$lat[idx]), method = "bilinear")
  # sst - -180-180 lon
  BGCArgo$sst[idx] = raster::extract(SST,cbind(BGCArgo$lon[idx], BGCArgo$lat[idx]), method = "bilinear")
  #ice - in a different projection
  si_points = data.frame(x = BGCArgo$lon[idx], y = BGCArgo$lat[idx])
  coordinates(si_points) <- ~x+y
  projection(si_points) <- "+proj=longlat +datum=WGS84"
  si_points = spTransform(si_points,SI@crs)
  BGCArgo$sea_ice[idx] = raster::extract(SI,si_points, method = "bilinear")
  
  # close files
  nc_close(fl1)
  # delte files
  unlink(file.path(dir, paste("CCMP_Wind_Analysis_",x,"_V02.0_L3.0_RSS.nc", sep = "")))
  # progress tracker
  print(which(u_date_list == x))
  # remove internal variables
  rm(idx, avWS, SST, file_location, file_location2)}

lon_ws[lon_ws<=0] = lon_ws[lon_ws<=0] + 360

write.table(BGCArgo,file.path(dirData,"BGCArgo60S_SO_sst_ice_31012020.tsv"),sep = "\t", row.names = F)
