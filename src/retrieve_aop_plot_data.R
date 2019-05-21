retrieve_aop_plot_data <- function(site, get_tile_chm = F){
  library(tidyverse)
  data <- readr::read_csv("./out/TOS_outputs/vegetation_structure_utm.csv") %>%
    filter(siteID %in% site)
  years <- readr::read_csv("./out/TOS_outputs/field_date_collection.csv") %>%
    filter(siteID %in% site)
  plots <- readr::read_csv("tmp/filesToStack10098/stackedFiles/vst_perplotperyear.csv") %>%
    filter(siteID %in% site) %>%
    select(plotID, easting, northing, utmZone) %>%
    unique
  
  
  #set paths to get data from
  source("./src/utilities.R")
  tileID <- unique(cbind(as.character(as.integer(plots$easting/1000)*1000),
                         as.character(as.integer(plots$northing/1000)*1000)))
  
  epsg <- get_epsg_from_utm(unique(years$siteID))
  
  paths <- get_aop_data_paths(fld = "/ufrc/ewhite/s.marconi/MMBRS/indir/AOP/",
                              NeonSites = site, year = years$scanDate, domainID= unique(data$domainID))
  if(get_tile_chm == T){
    crownITC(paths$pt, wd = "./AOP_from_coords/", 
             pttrn = paste(tileID[,1], "_", tileID[,2], sep=""),
             epsg = epsg,  chm_f = paths$chm_f, dtm_pt = paths$dtm_pt,
             pybin = "/home/s.marconi/.conda/envs/quetzal3/bin")
  }
  #get canopy height from lidar
  aop_chm_plot(plots, tileID, epsg, paths)
  #get itcs from lidar
  aop_itcs_plot()
  #get hiperspectral
  aop_hiperspectral()
}