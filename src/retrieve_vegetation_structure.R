retrieve_vegetation_structure <- function(){
  #'
  #'
  #'
  file_mapping = read_csv("./tmp/filesToStack10098/stackedFiles/vst_mappingandtagging.csv") %>%
    dplyr::select(c("individualID", "eventID", "domainID","siteID","plotID","subplotID",
                    "nestedSubplotID","pointID","stemDistance","stemAzimuth",
                    "supportingStemIndividualID","previouslyTaggedAs",
                    "taxonID","scientificName"))

  plots<-sf::st_read("./dat/TOS_inputs/All_Neon_TOS_Points_V5.shp") %>% filter(str_detect(appMods,"vst"))
  
  
  
  dat<-file_mapping %>% 
    mutate(pointID=factor(pointID, levels = levels(unique(plots$pointID))) )%>% 
    mutate(plotID=factor(plotID, levels = levels(unique(plots$plotID)))) %>% 
    inner_join(plots,by=c("plotID","pointID"))
  
  dat <- dat[!is.na(dat$stemAzimuth), ]
  # dat <- dat %>% filter(siteID %in% c("DSNY", "GRSM", "GUAN", "HARV",  "KONZ" ,
  #                                     "LENO", "MLBS", "MOAB", "ORNL", "SCBI", 
  #                                     "SERC",  "STEI", "TALL", "UKFS"))
  # # get tree coordinates
  dat_apply <- dat %>%
    dplyr::select(c(stemDistance, stemAzimuth, easting, northing)) 
  coords <- apply(dat_apply,1,function(params)retrieve_dist_to_utm(params[1],params[2], params[3], params[4])) %>%
    t %>%
    data.frame
  colnames(coords) <- c('UTM_E', 'UTM_N')
  field_tag <- cbind(dat, coords) %>% filter(!is.na(UTM_E))
  
  max_no_na <- function(x)max(x, na.rm=T)
  apparent = readr::read_csv("./tmp/filesToStack10098/stackedFiles/vst_apparentindividual.csv") %>%
    dplyr::select("individualID", "stemDiameter", "height", "maxCrownDiameter",  
                  "ninetyCrownDiameter", "growthForm", "plantStatus", "canopyPosition", "shape") %>%
    group_by(individualID) %>%
    summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.)))
  
  apparent$stemDiameter[is.infinite(apparent$stemDiameter)] <- NA
  apparent$height[is.infinite(apparent$height)] <- NA
  apparent$maxCrownDiameter[is.infinite(apparent$maxCrownDiameter)] <- NA
  apparent$ninetyCrownDiameter[is.infinite(apparent$ninetyCrownDiameter)] <- NA
  
  crown_attributes = left_join(field_tag, apparent, by="individualID") %>%
    unique
  
  #check for dopplegangers
  which_multiple <- crown_attributes %>% group_by(individualID) %>% 
    summarize(n=n()) %>% filter(n>1)  
  
  #find a more elegant solution for dopplegangers (at the time of Aug 2019 mainly just NEON.PLA.D05.STEI.01071A and NEON.PLA.D14.SRER.01112)
  # crown_attributes %>% group_by(individualID) %>% 
  # top_n(1, wt = individualID)
  
  crown_attributes <- get_crown_dimensions(crown_attributes)
  write_csv(crown_attributes, './out/TOS_outputs/field_data.csv')
}
