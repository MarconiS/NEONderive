retrieve_aop_plot_data <- function(site){
  data <- readr::read_csv("./out/TOS_outputs/vegetation_structure_utm.csv") %>%
    filter(siteID %in% site)
  years <- read_csv("./out/TOS_outputs/field_date_collection.csv") %>%
    filter(siteID %in% site)
  plots <- read_csv("tmp/filesToStack10098/stackedFiles/vst_perplotperyear.csv") %>%
    filter(siteID %in% site) %>%
    select(plotID, easting, northing, utmZone) %>%
    unique
  
  
  #set paths to get data from
  tileID <- unique(cbind(as.character(as.integer(plots$easting/1000)*1000),
                         as.character(as.integer(plots$northing/1000)*1000)))
  
  epsg <- get_epsg_from_utm(unique(centroids$siteID))
  
  paths <- get_aop_data_paths(fld = "/ufrc/ewhite/s.marconi/MMBRS/indir/AOP/DP1.30003.001/")
  #get canopy height from lidar
  aop_chm_plot()
  #get itcs from lidar
  aop_itcs_plot()
  #get hiperspectral
  aop_hiperspectral()
}