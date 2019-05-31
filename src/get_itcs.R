library(lidR)
library(tidyverse)
library(raster)
library(rgdal)  # input/output, projections
library(sf)
chm_rgba = raster("~/Documents/Data/plot/false_col/MLBS_072_out.tif")
chm = raster("~/Documents/Data/plot/CHM/CHMMLBS_072.tif")
las = 
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

rescale <- function(x, max_val) (x-min(as.matrix(x), na.rm=T))/
  (max(as.matrix(x), na.rm=T) - min(as.matrix(x), na.rm=T)) * max_val
chm_rgba = rescale(x = chm_rgba, max_val = 50)

# treetops <- data %>% dplyr::filter(plotID == "MLBS_072") %>%
#   dplyr::select("individualID", "UTM_E", "UTM_N", "height", "stemDiameter") %>% #, 
#   unique %>%
#   dplyr::group_by(UTM_E, UTM_N) %>%
#   summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.))) %>%
#   st_as_sf(coords = c("UTM_E", "UTM_N"), crs = epsg)

source("./src/utilities.R")
settbuff=buffer_bbox(st_bbox(treetops),3)
attr(settbuff, "class") = "bbox"
attr(st_geometry(treetops), "bbox") = settbuff
#crop the chm to the extent of field data collection 
#chm_itc <- raster::crop(chm_rgb, extent(treetops))
chma_itc <- raster::crop(chm_rgba, extent(treetops))
summary(chma_itc)
#use tree data allometry to apply dalponte 2016
#itcs
unsupervised = tree_detection(chm_itc,  lmf(ws = 5))
ttops <- tree_detection(chm_rgba, lmf(2,2))
plot(chm_rgba)
plot(ttops, add = T)
#algo = dalponte2016(chm_itc, treetops = unsupervised)()

# dalpo = dalponte2016(chm_rgba, treetops = as(treetops, "Spatial"), 
#                      th_cr = 0.4, th_seed = 0.4)()
# water = mcwatershed(chm_rgba, treetops = as(treetops, "Spatial"))()
# #mcwatershed(chm, treetops, th_tree = 2, ID = "treeID")
# silva = silva2016(chm_rgba, treetops = as(treetops, "Spatial"))()

dalpo = dalponte2016(chm_rgba, treetops = ttops)()
water = mcwatershed(chm_rgba, treetops = ttops)()
#mcwatershed(chm, treetops, th_tree = 2, ID = "treeID")
silva = silva2016(chm_rgba, treetops = ttops)()

plot(dalpo, col=sample(color, max(as.matrix(dalpo), na.rm=T),))
plot(water, col=sample(color, max(as.matrix(dalpo), na.rm=T)))
plot(silva, col=sample(color, max(as.matrix(dalpo), na.rm=T)))
crs(dalpo) = crs(water) = crs(silva) = crs(chm_itc)

algo = dalponte2016(chm_itc, treetops = ttops)
las_dp  = lastrees(las, algo)
metric = tree_metrics(las_dp, .stdtreemetrics)
hulls  = tree_hulls(las_dp, attribute = "treeID")
hulls@data = dplyr::left_join(hulls@data, metric@data)

proj4string(hulls) <- crs(treetops)
itcs_dp <- st_join(st_as_sf(hulls), treetops) 
st_write(itcs_dp, paste("./out/AOP/plot/ITCs/", plots[jj,"plotID"], ".shp"), delete_layer=TRUE)

algo = silva2016(chm_itc, treetops = ttops)
las_sl  = lastrees(las, algo)
metric = tree_metrics(las_sl, .stdtreemetrics)
hulls  = tree_hulls(las_sl, attribute = "treeID")
hulls@data = dplyr::left_join(hulls@data, metric@data)

proj4string(hulls) <- crs(treetops)
itcs_sl <- st_join(st_as_sf(hulls), treetops)
st_write(itcs_sl, paste("./out/AOP/plot/ITCs/", plots[jj,"plotID"], ".shp"), delete_layer=TRUE)

algo = dalponte2016(chm_itc, treetops = ttops)
las_ws  = lastrees(las, algo)
metric = tree_metrics(las_ws, .stdtreemetrics)
hulls  = tree_hulls(las_ws, attribute = "treeID")
hulls@data = dplyr::left_join(hulls@data, metric@data)

#add itc names and save itcs
proj4string(hulls) <- crs(treetops)
itcs_ws <- st_join(st_as_sf(hulls), treetops)
st_write(itcs_ws, paste("./out/AOP/plot/ITCs/", plots[jj,"plotID"], ".shp"), delete_layer=TRUE)



# crs(dalpo) = crs(chm_itc)
# plot(chm_rgba)
# plot(st_geometry(treetops), add = T)
# plot(ttops, add = T)
# plot(chm_rgba) 
# plot(dalpo,  col=sample(color, max(as.matrix(dalpo), na.rm=T)))
# plot(st_geometry(treetops), add = T, fill = 'black')
