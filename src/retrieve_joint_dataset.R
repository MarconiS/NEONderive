retrieve_joint_dataset <- function(){
  #' clean leaf and structural traits data and create final dataset, 
  #' utm_dataset.csv including leaf traits,
  #' field_traits_dataset.csv including only data points for which both 
  #' coordinates, vegetation structure, and leaf traits  available
  #' 
  #'@import dplyr

  chemical <- read_csv("./out/TOS_outputs/chemical_data.csv") %>%
    select(-one_of("toxicodendronPossible.x", "toxicodendronPossible.y", 
                   "chlorophyllSampleCode.x", "chlorophyllSampleCode.y",
                   "ligninSampleBarcode.x", "ligninSampleBarcode.y", 
                   "cnSampleCode.x", "cnSampleCode.y"))
  
  date <- chemical %>% select(scanDate, siteID) %>% unique()
  date$scanDate<- lubridate::year(date$scanDate)
  write_csv(unique(date), "./out/TOS_outputs/field_date_collection.csv")
  
  # isotope <- readr::read_csv("./TOS_retriever/out/isotopes_data.csv") %>%
  #   select(sampleID, d13C, d15N)
  structure <- read_csv("./out/TOS_outputs/field_data.csv") %>%
    select(-one_of("plotType", "coordinateUncertainty",
                   "subplotID"))
  
  #uniform TOS with AOP reference systems. Converting STEI into STEI and CHEQ
  source("./src/utilities.R")
  new_cords_stei <- convert_stei(structure[which(structure$siteID=="STEI"),])
  structure[which(structure$siteID=="STEI"),c("siteID", "UTM_E", "UTM_N")] <- 
    new_cords_stei %>% select(siteID, easting, northing)
  write_csv(structure, "./out/TOS_outputs/vegetation_structure_utm.csv")
  structure <- structure %>% select(-one_of("elevation"))
  
  ft_nm <- c("leafMassPerArea", "dryMassFraction", "ligninPercent", "cellulosePercent", "foliarPhosphorusConc","foliarPotassiumConc", "foliarCalciumConc",  "extractChlAConc","extractChlBConc","extractCarotConc",  
             "foliarMagnesiumConc","foliarSulfurConc","foliarManganeseConc","foliarIronConc","foliarCopperConc",
             "foliarBoronConc","foliarZincConc", "nitrogenPercent","carbonPercent")
  
  chemical <- chemical[!is.na(chemical$individualID),]
  chemical[ft_nm] <- apply(chemical[ft_nm], 2,  function(x)remove_outliers(x))
  chemical <- chemical %>% group_by(.dots=c("individualID")) %>% 
    summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.))) %>% 
    filter("subsample1Height" > 3) %>%
    select("individualID","domainID", "taxonID", "subsample1Height","scientificName", 
           "elevation", "siteID", "decimalLatitude","decimalLongitude", "leafMassPerArea", 
           "ligninPercent", "cellulosePercent", "foliarPhosphorusConc",
           "foliarPotassiumConc", "foliarCalciumConc", "extractChlAConc","extractChlBConc",
           "extractCarotConc",  "foliarMagnesiumConc","foliarSulfurConc",
           "foliarManganeseConc","foliarIronConc","foliarCopperConc",
           "foliarBoronConc","foliarZincConc", "nitrogenPercent","carbonPercent") 
  
  ggplot(tidyr::gather(chemical[-c(1:3,5:9)]), aes(value)) + 
    geom_histogram(bins = 40) + geom_rug()+
    facet_wrap(~key, scales = 'free')
  
  chemical = chemical[complete.cases(chemical),] %>% unique
  readr::write_csv(chemical, './out/TOS_outputs/tree_traits_dataset.csv')
  
  #chemical <- chemical %>% select(-one_of("siteID"))
  # just the geolocalized data
  full_data <- inner_join(chemical, structure, 
      by = c("individualID","siteID","domainID", "taxonID", "scientificName")) %>% 
    unique %>%
    group_by(individualID) %>%
    summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.)))

  readr::write_csv(full_data, './out/TOS_outputs/utm_dataset.csv')
  
}
