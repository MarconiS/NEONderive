
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.05, .95), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

get_mod_r2 <- function(pred, obs){
  #1 - sum((pred - obs)^2) / sum((obs - mean(obs, na.rm=T))^2)
  1 - sum((pred - obs)^2) / sum((obs - mean(obs, na.rm=T))^2)
} 

# OutVals = apply(y_obs[-1], 2,  function(x)remove_outliers(x))
# summary(OutVals)


plot_spectra<-function(plt_dat){
  #plot reflectances
  plot_data <- plt_dat %>% 
    dplyr::select(-one_of(c("siteID", "taxonID",  "band_site","band_species", "flightpath"))) 
  plot_data <- plot_data %>% dplyr::select(-one_of("individualID")) %>%
    t %>%
    data.frame
  colnames(plot_data) = unlist(plt_dat[1]) # the first row will be the header
  plot_data <- data.frame(bnd = 1:dim(plot_data)[1], plot_data)
  ggdat <- tidyr::gather(plot_data, treeID, Reflectance,-bnd)
  
  return(ggplot(ggdat, aes(x = bnd, y = Reflectance)) + 
           geom_line(aes(color = factor(treeID), alpha= 1), size = 0.2) +
           theme_bw()+
           theme(legend.position="none"))
  
}

convert_stei <- function(dat){
  library(rgdal)
  what_to_keep <- which(dat$easting > 270000)
  tmp_keep <- dat[what_to_keep,]
  dat$utmZone <- "15N"
  dat$siteID <- "CHEQ"
  coordinates(dat) <- c("easting", "northing")
  proj4string(dat) <- CRS("+init=epsg:32616") # WGS 84
  dat <- spTransform(dat, CRS("+init=epsg:32615"))
  new_dat <- cbind(dat@data, dat@coords)
  new_dat[what_to_keep, "easting"] <- new_dat[what_to_keep, "UTM_E"]
  new_dat[what_to_keep, "northing"] <- new_dat[what_to_keep, "UTM_N"]
  new_dat[what_to_keep, "siteID"] <- tmp_keep[, "siteID"]
  
  return(new_dat)
}

get_aop_data_paths <- function(fld, year, domainID, NeonSites){
  
  if(NeonSites %in% c("SERC", "GRSM", "SCBI", "CHEQ", "STEI", "TALL")){
    year = 2017
  }
  pt <- paste(fld, "/DP1.30003.001/", year, "/FullSite/", unique(domainID), "/",
              "/", year,"_", NeonSites, "_2/", "L1/DiscreteLidar/ClassifiedPointCloud/", sep="")
  f_path <- paste(fld, "/DP3.30006.001/", year, "/FullSite/", unique(domainID), "/",
                  "/", year,"_", NeonSites,"_2/", "L3/Spectrometer/Reflectance/", sep="") #H5
  chm_f <- paste(fld, "/DP1.30003.001/", year, "/FullSite/", unique(domainID), "/",
                 "/", year,"_", NeonSites, "_2/", "L3/CHM/", sep="")
  dtm_pt = paste(fld, "/DP3.30024.001/", year, "/FullSite/", unique(domainID), "/",
                 "/", year,"_", NeonSites, "_2/", "L3/DiscreteLidar/DTMGtif/", sep="")
  
  rgb_pt = paste(fld, "/DP3.30010.001/", year, "/FullSite/", unique(domainID), "/",
                 "/", year,"_", NeonSites, "_2/", "L3/DiscreteLidar/DTMGtif/", sep="")
  
  return(list(pt = pt, f_path = f_path, chm_f = chm_f, dtm_pt= dtm_pt))
  
}

get_downloaded_data_paths <- function(fld = "//orange/ewhite/NeonData/", year, domainID, NeonSites){
  
  if(NeonSites %in% c("SERC", "GRSM", "SCBI", "CHEQ", "STEI", "TALL")){
    year = 2017
  }
  recurred = list.files(paste(fld, NeonSites,"/DP1.30003.001/", year, "/FullSite/", unique(domainID), "/", sep=""))
  recurred = substr(recurred, nchar(recurred), nchar(recurred))
  pt <- paste(fld, NeonSites,"/DP1.30003.001/", year, "/FullSite/", unique(domainID), "/",
              "/", year,"_", NeonSites, "_", recurred, "/", "L1/DiscreteLidar/ClassifiedPointCloud/", sep="")
  f_path <- paste(fld, NeonSites, "/DP3.30006.001/", year, "/FullSite/", unique(domainID), "/",
                  "/", year,"_", NeonSites,"_", recurred, "/", "L3/Spectrometer/Reflectance/", sep="") #H5
  chm_f <- paste(fld,NeonSites, "/DP1.30003.001/", year, "/FullSite/", unique(domainID), "/",
                 "/", year,"_", NeonSites, "_", recurred, "/", "L3/CHM/", sep="")
  dtm_pt = paste(fld,NeonSites, "/DP3.30024.001/", year, "/FullSite/", unique(domainID), "/",
                 "/", year,"_", NeonSites, "_", recurred, "/", "L3/DiscreteLidar/DTMGtif/", sep="")
  
  rgb_pt = paste(fld, NeonSites,"/DP3.30010.001/", year, "/FullSite/", unique(domainID), "/",
                 "/", year,"_", NeonSites, "_", recurred, "/", "L3/Camera/Mosaic/V01/", sep="")
  
  return(list(pt = pt, f_path = f_path, chm_f = chm_f, dtm_pt= dtm_pt, rgb_pt = rgb_pt))
  
}
# get_epsg_from_utm <- function(utm){
#   dictionary <- cbind(32616, 32618, 32615, 32617, 32616, 32616, 32616, 32612, 32613, 32617, 32617, 
#                       32614, 32618, 32616, 32619, 32617, 32615) 
#   colnames(dictionary) <- c("STEI", "SERC","CHEQ", "SCBI", "GRSM", "ORNL", "TALL", "MOAB", 
#                             "JORN", "OSBS", "MLBS", "KONZ", "HARV", "LENO", "GUAN", "DSNY", "UKFS")
#   return(dictionary[colnames(dictionary)==utm])
# }

buffer_bbox <- function(bbox, buffer){
  bbox[c(1,2)] <- bbox[c(1,2)] - buffer
  bbox[c(3,4)] <- bbox[c(3,4)] + buffer
  return(bbox)
}

get_itcs_in_tile <- function(hps_fi = NULL, centroids = NULL, f_path = NULL, NeonSites = NULL, buffer = 2){
  library(rhdf5)
  h5name <- paste(f_path,hps_fi,sep = "/")
  #mapInfo <- h5read(h5name,"map info")
  #mapInfo<-unlist(strsplit(mapInfo, ","))
  
  #new metadata ['Metadata']['Coordinate_System']['Map_Info']
  mapInfo <- h5read(h5name,paste(NeonSites, "/Reflectance/Metadata", sep=""))
  mapInfo <-unlist(strsplit(mapInfo$Coordinate_System$Map_Info, ","))
  #grab the utm coordinates of the lower left corner
  xMin<-as.numeric(mapInfo[4])
  yMax<-as.numeric(mapInfo[5])
  reflInfo <- h5ls(h5name)
  #reflInfo <- unlist(strsplit(reflInfo$dim[3], split = "x"))
  reflInfo <- unlist(strsplit(reflInfo$dim[7], split = "x"))
  
  nRows <- as.integer(reflInfo[2])
  nCols <- as.integer(reflInfo[1])
  xMax <- xMin + nCols
  yMin <- yMax - nRows
  #be sure you canclip a 60 by 60 pic
  mask.x <- (centroids$easting < xMax-buffer) & (centroids$easting > buffer + xMin)
  mask.y <- (centroids$northing < yMax-buffer) & (centroids$northing > yMin+buffer)
  itcextract <- centroids[mask.x & mask.y,]
  #h5closeAll()
  return(itcextract)
}

utm_to_geographic <- function(dat){
  # from UTM zone and E N coords, get the geographic lat and longitude
  latlon = NULL
  library(sf)
  for(ii in unique(dat$utmZone)){
    epsg = paste("+init=epsg:326", gsub('.{1}$', '', ii), sep="")
    foo = dat %>% filter(utmZone == ii)
    foo = sf::st_as_sf(foo, coords = c("UTM_E", "UTM_N"), crs = epsg)
    foo = st_transform(foo, "+init=epsg:4326")
    
    latlon = rbind.data.frame(latlon, cbind(foo$unique_id, st_coordinates(foo)),  stringsAsFactors = F)
  }
  colnames(latlon) <- c("unique_id", "longitude", "latitude")
  new_pneon <- inner_join(dat, data.frame(latlon), by ="unique_id")
  write_csv(new_pneon, "./out/preliminary_neon_points.csv")
  return(latlon)
  
}

get_epsg_from_utm <- function(utm){
  utm <-  substr(utm,1,nchar(utm)-1)
  if(as.numeric(utm)<10)utm <- paste('0', utm, sep="")
  epsg <- paste("326", utm, sep="")
  return(epsg)
}

get_lat_long <- function(field_data){
  library(rgdal)
  new_dat <- NULL
  for(NeonSites in unique(field_data$siteID)){
    dat <- field_data[field_data$siteID %in% NeonSites,] 
    epsg <- get_epsg_from_utm(unique(dat$utmZone))
    dat <- dat[complete.cases(dat$UTM_E), ]
    utm_coords <- dat[c("UTM_E", "UTM_N")]
    coordinates(dat) <- c("UTM_E", "UTM_N")
    proj4string(dat) <- CRS(paste("+init=epsg:", epsg, sep ="")) 
    CRS.new <- CRS("+init=epsg:4326")
    dat <- spTransform(dat, CRS.new)
    coords_dat <- dat@coords
    colnames(dat@coords) <- c("itc_lat", "itc_lon")
    new_dat <- rbind(new_dat, cbind(dat@data, utm_coords, dat@coords))
  }
  return(new_dat)
}
