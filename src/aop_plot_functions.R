aop_chm_plot <- function(plots, tileID, epsg, paths, bff = 25, cl = 16){
  #' use lidR to clip the lidar data, create a very high resolution chm, and return predicted itcs
  #' given we know the tree tops from NEON TOS
  #' @include lidR
  #' @param plots data.frame.
  #' 
  #' @param tileID data.frame.
  #' @param epsg 
  #' @param paths   
  
  library(lidR)
  library(sf)
  dir.create(file.path("out", "AOP"))#, showWarnings = FALSE)
  dir.create(file.path("out", "AOP", "plot"))#, showWarnings = FALSE)
  dir.create(file.path("out", "AOP", "plot", "CHM"))#, showWarnings = FALSE)
  dir.create(file.path("out", "AOP", "plot", "ITCs"))#, showWarnings = FALSE)
  dir.create(file.path("out", "AOP", "spectra"))#, showWarnings = FALSE)
  
  library(foreach)
  library(doParallel)
  registerDoSEQ()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  clusterCall(cl, function(x) .libPaths(x), .libPaths())
  
  results <- foreach(ii = 1:nrow(tileID)) %dopar% {
  #for(ii in 1:nrow(tileID)){
    lid_tile <- list.files(paths$pt, pattern = paste(tileID[ii,], collapse = "_"))
    las_tl <- readLAS(paste(paths$pt, lid_tile, sep="/"))
    hps_f = list.files("//ufrc/ewhite/s.marconi/MMBRS/AOP_from_coords/outputs/itcTiff", pattern = paste(tileID[ii,], collapse = "_"))
    hps <- raster::stack(paste(paths$f_path, hps_f, sep="/"))
    
    for(jj in 1:nrow(plots)){
      crds <- plots[jj, c("easting", "northing")]
      las <- lasclipRectangle(las_tl, xleft = crds$easting - bff, 
                              crds$northing - bff, 
                              crds$easting + bff, 
                              crds$northing + bff)
      laspl = lasnormalize(las, tin())
      thr <- c(0,2,5,10,15)
      edg <- c(0, 1.5)
      chm <- grid_canopy(laspl, 0.25, pitfree(thr, edg))
      chm <- stretch(chm, minq=0.05, maxq=0.95)
      
      raster::writeRaster(chm, filename=paste("./out/AOP/plot/CHM",
                                              plots[jj,"plotID"], ".tif", sep=""), 
                          format="GTiff", overwrite=TRUE)
      treetops <- data %>% filter(plotID == unlist(plots[jj,"plotID"])) %>%
        dplyr::select("individualID", "UTM_E", "UTM_N", "height", "stemDiameter") %>% #, 
        #       "maxCrownDiameter", "ninetyCrownDiameter") %>%
        unique %>%
        group_by(UTM_E, UTM_N) %>%
        summarise_all(list(~if(is.numeric(.)) mean(., na.rm = TRUE) else first(.))) %>%
        st_as_sf(coords = c("UTM_E", "UTM_N"), crs = epsg)
      
      #use plot extent plus 3m buffer
      settbuff=buffer_bbox(st_bbox(treetops),3)
      attr(settbuff, "class") = "bbox"
      attr(st_geometry(treetops), "bbox") = settbuff
      #crop the chm to the extent of field data collection 
      chm_itc <- raster::crop(chm, extent(treetops))
      #use tree data allometry to apply dalponte 2016
      #itcs
      algo = dalponte2016(chm_itc, treetops = as(treetops, "Spatial"))
      las  = lastrees(las, algo)
      metric = tree_metrics(las, .stdtreemetrics)
      hulls  = tree_hulls(las)
      hulls@data = dplyr::left_join(hulls@data, metric@data)
      
      #add itc names and save itcs
      proj4string(hulls) <- crs(treetops)
      itcs <- st_join(st_as_sf(hulls), treetops) %>%
        filter(!is.na(individualID))
      st_write(itcs, paste("./out/AOP/plot/ITCs/", plots[jj,"plotID"], ".shp"))
      
      #extract data from hiperspectral
      spectra <- extract(hps, itcs, df=TRUE) 
      write_csv(spectra, paste("./out/AOP/plot/spectra/", plots[jj,"plotID"], ".csv"))
      
    }
  }
}
