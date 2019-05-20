
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
  
  pt <- paste(fld, "/DP1.30003.001/", 
              year, "/FullSite/", unique(centroids$domainID), "/",
              "/", year,"_", NeonSites, "_2/", "L1/DiscreteLidar/ClassifiedPointCloud/", sep="")
  f_path <- paste(fld, "/DP3.30006.001/", year, "/FullSite/", unique(domainID), "/",
                  "/", year,"_", NeonSites,"_2/", "L3/Spectrometer/Reflectance/", sep="") #H5
  chm_f <- paste(fld, "/DP1.30003.001/", year, "/FullSite/", unique(domainID), "/",
                 "/", year,"_", NeonSites, "_2/", "L3/CHM/", sep="")
  dtm_pt = paste(fld, "/AOP/DP3.30024.001/", year, "/FullSite/", unique(domainID), "/",
                 "/", year,"_", NeonSites, "_2/", "L3/DiscreteLidar/DTMGtif/", sep="")
  
  return(list(f_path, chm_f, dtm_pt))
  
}
