rs_data_main <- function(site = NULL, get_tile_chm = T, getAOP = F){
  library(tidyverse)
  #years = data.frame(scanDate = 2018, siteID = "OSBS", stringsAsFactors = F)
  years <- readr::read_csv("./out/TOS_outputs/field_date_collection.csv") %>%
    dplyr::filter(siteID %in% site)
  
  data <- readr::read_csv("./out/TOS_outputs/vegetation_structure_utm.csv") %>%
    dplyr::filter(siteID %in% site)
  plots <- readr::read_csv("tmp/filesToStack10098/stackedFiles/vst_perplotperyear.csv") %>%
    dplyr::filter(siteID %in% site) %>%
    dplyr::select(plotID, easting, northing, siteID, utmZone) %>%
    unique
  
  #set paths to get data from
  source("./src/utilities.R")
  file.sources = paste("./src/",
                       list.files("./src/", pattern="aop"), sep="/") 
  sapply(file.sources,source,.GlobalEnv)
  
  tileID <- unique(cbind(as.character(as.integer(plots$easting/1000)*1000),
                         as.character(as.integer(plots$northing/1000)*1000)))
  
  epsg <- get_epsg_from_utm(unique(years$siteID))
  
  paths <- get_aop_data_paths(fld = "/ufrc/ewhite/s.marconi/MMBRS/indir/AOP/",
                              NeonSites = site, year = years$scanDate, domainID= unique(data$domainID))
  if(get_tile_chm == T){
    dir.create(file.path("out", "AOP", "tiles"))#, showWarnings = FALSE)
    dir.create(file.path("out", "AOP", "plot", "itcTiff"))#, showWarnings = FALSE)
    
    #get AOP data, if needed
    if(getAOP){
      aop_retrieve(fin = plots)
    }
    
    centroids <- plots
    centroids$taxonID <- "NA"
    centroids$individualID <- centroids$plotID
    #
    aop_canopy_height(pt = paths$pt, wd = "./",
             pttrn = paste(tileID[,1], "_", tileID[,2], sep=""),
             epsg = epsg,  chm_f = paths$chm_f, dtm_pt = paths$dtm_pt, cores = 4,
             pybin = "/home/s.marconi/.conda/envs/quetzal3/bin")

    hps_f = list.files(paths$f_path)
    aop_hps_data(centroids = centroids, hps_f = hps_f, f_path =  paths$f_path, buffer = 25,
                       chm_f = paths$chm_f, epsg=epsg, wd =  "./", NeonSites=site, cores = 5)
    
  }
  
  #get canopy height from lidar
  aop_chm_plot(plots, data, tileID, epsg, paths, bff = 25, cores = 4)

}