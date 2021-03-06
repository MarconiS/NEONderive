#get neon TOS
aop_retrieve <- function(fin){
  library(neonUtilities)
  library(tidyverse)
  library(readr)
  #fin <- readr::read_csv("./TOS_retriever/out/field_data_utm.csv")
  years <- read_csv("./out/TOS_outputs/field_date_collection.csv") %>%
    dplyr::filter(siteID == unique(fin$siteID))
  tryCatch({
  years[years$siteID == "GRSM", "scanDate"] <- 2017
  years[years$siteID == "SERC", "scanDate"] <- 2017
  })
  #tiles <- coords_for_tiles[-1] %>% unique %>% filter((siteID %in% c("ORNL")))
  
  #coords_for_tiles <- fin %>% dplyr::select(individualID, siteID, utmZone, UTM_E, UTM_N) %>%
  coords_for_tiles <- fin %>% dplyr::select(plotID, siteID, utmZone, easting, northing) %>%
    dplyr::filter(siteID %in% years$siteID)
  coords_for_tiles$easting <- as.integer(coords_for_tiles$easting / 1000) * 1000
  coords_for_tiles$northing <- as.integer(coords_for_tiles$northing / 1000) * 1000
  
  tiles <- coords_for_tiles[-1] %>% unique %>% filter(!(siteID %in% c("SOAP", "DEJU", "ORNL")))
  tiles[tiles$siteID == "CHEQ", "siteID"] <- "STEI"
  tiles %>% unique
  #elevation
  for(ii in 1:nrow(tiles)){
    
    #elevation
    byTileAOP("DP3.30024.001", site = tiles[ii,"siteID"],
              year = years[years$siteID %in% tiles[ii, "siteID"], "scanDate"],
              tiles[ii,"easting"], tiles[ii,"northing"],
              buffer = 0, check.size = F, savepath = "./indir/AOP/")
    # #reflectance
    byTileAOP("DP3.30006.001", site = tiles[ii,"siteID"],
              year = years[years$siteID %in% tiles[ii, "siteID"], "scanDate"],
              tiles[ii,"easting"], tiles[ii,"northing"],
              buffer = 0, check.size = F, savepath = "./indir/AOP/")
    #aspect and slope
    byTileAOP("DP3.30025.001", site = tiles[ii,"siteID"],
              year = years[years$siteID %in% tiles[ii, "siteID"], "scanDate"],
              tiles[ii,"easting"], tiles[ii,"northing"],
              buffer = 0, check.size = F, savepath = "./indir/AOP/")

    #rgb ortho
    byTileAOP("DP3.30010.001", site = tiles[ii,"siteID"], 
              year = years[years$siteID %in% tiles[ii, "siteID"], "scanDate"], 
              tiles[ii,"easting"], tiles[ii,"northing"],
              buffer = 0, check.size = F, savepath = "./indir/AOP/")
    #1m2 chm
    byTileAOP("DP3.30015.001", site = tiles[ii,"siteID"], 
              year = years[years$siteID %in% tiles[ii, "siteID"], "scanDate"], 
              tiles[ii,"easting"], tiles[ii,"northing"],
              buffer = 0, check.size = F, savepath = "./indir/AOP/")
    
    
  }
  
  # #L1_2_dat <- years %>% filter((siteID %in% unique(tiles$siteID)))
  # byFileAOP("DP1.30003.001", site = years$siteID, year = years$scanDate, 
  #             check.size = F, savepath = "./indir/AOP/")
}


