## This code matches SOCAT wind speed data

library(data.table)
library(RCurl)
library(ggplot2)
library(ncdf4)
library(cwhmisc)
library(stringr)

dirData = "/perm_storage/home/kbaldry/Lennart/Data"
SOCAT = fread(file.path(dirData, "SOCATv2_SouthernOceans.tsv"), skip = "/starthere", header = T, sep = "\t")
colnames(SOCAT)[10:11] = c("lon","lat")
SOCAT$lon[SOCAT$lon>180] = SOCAT$lon[SOCAT$lon>180] - 360
south_SOCAT = SOCAT %>% filter(lat < -60, fCO2rec_flag < 4)
write.table(south_SOCAT,file.path(dirData,"SOCAT_subset.tsv"), row.names = F,sep = "\t")
# south_SOCAT = fread(file.path(dirData, "SOCATv2_SO_ws_dailyaverage.tsv"), header = T, sep = "\t")


url = "http://data.remss.com/ccmp/v02.0"
url2 = "https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/access/avhrr-only/"

south_SOCAT$date = as.POSIXct(paste(south_SOCAT$yr,sprintf("%02.0f",south_SOCAT$mon), sprintf("%02.0f",south_SOCAT$day),sprintf("%02.0f",south_SOCAT$hh),sprintf("%02.0f",south_SOCAT$mm),sprintf("%02.0f",south_SOCAT$ss) , sep = ""), format = "%Y%m%d%H%M%S")


date_list = paste(south_SOCAT$yr,sprintf("%02.0f",south_SOCAT$mon), sprintf("%02.0f",south_SOCAT$day),sep = "")
u_date_list = unique(date_list)
dir = "/perm_storage/home/kbaldry/Lennart/Data/Wind_data"
dir2 = "/perm_storage/home/kbaldry/Lennart/Data/SST_data"
u_date_list = u_date_list[-which(as.numeric(u_date_list)<19870000)]
south_SOCAT$ws = NA
south_SOCAT$sst = NA
south_SOCAT$sea_ice = NA
# u_date_list2 = u_date_list[c(seq(1,901,8),seq(2,901,8),seq(3,901,8),seq(4,901,8),seq(5,901,8),seq(6,901,8),seq(7,901,8),seq(8,901,8))]
# # if the length of ids is 0 after this, then stop
# numJobs <- length(u_date_list2)
# if (numJobs > 0) {
#   
#   numCores <- min(numJobs, 8)
#   
#   theCluster <- parallel::makeForkCluster(getOption("cl.cores", numCores))
#   
#   discard_output <- parallel::clusterApplyLB(theCluster
#                                                ,x=u_date_list2
#                                                ,fun=Get_ws)
#   
#   
#   parallel::stopCluster(theCluster)
#   
# }
# 



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
  
  # find index of mathing dates
  idx = which(date_list == x)
  # Bilinear interpolation
  # wind speed - 0-360 lon
  lon_ws = south_SOCAT$lon[idx]
  lon_ws[lon_ws<=0] = lon_ws[lon_ws<=0] + 360
  south_SOCAT$ws[idx] = raster::extract(avWS,cbind(lon_ws, south_SOCAT$lat[idx]), method = "bilinear")
  # sst - -180-180 lon
  south_SOCAT$sst[idx] = raster::extract(SST,cbind(south_SOCAT$lon[idx], south_SOCAT$lat[idx]), method = "bilinear")
  #ice - in a different projection
  si_points = data.frame(x = BGCArgo$lon[idx], y = BGCArgo$lat[idx])
  coordinates(si_points) <- ~x+y
  projection(si_points) <- "+proj=longlat +datum=WGS84"
  si_points = spTransform(si_points,SI@crs)
  south_SOCAT$sea_ice[idx] = raster::extract(SI,si_points, method = "bilinear") 
  
  # close files
  nc_close(fl1)
  # delte files
  unlink(file.path(dir, paste("CCMP_Wind_Analysis_",x,"_V02.0_L3.0_RSS.nc", sep = "")))
  # progress tracker
  print(which(u_date_list == x))
  # remove internal variables
  rm(idx, avWS, SST,SI, file_location, file_location2)}

write.table(south_SOCAT,file.path(dirData,"SOCATv2_SO_sst_ice_31012020.tsv"),sep = "\t", row.names = F)