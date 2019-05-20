retrieve_dist_to_utm <- function(dista, angle, xcord, ycord){
  #' Using the coordinates of the point ids and trigonometry, 
  #' get precise locations of trees
  #'@param dista numeric. distance between point id and tree stem
  #'@param angle numeric. Azimut angle between point id and tree stem
  #'@param xcord numeric. UTM East of the point ID
  #'@param ycord numeric. UTM West of the point ID
  #'
  #' @export point.coords numeric. UTM E and N of the individual tree
  #' @examples 
  #' if (interactive()) {
  #'   addinslistAddin()
  #' }
  #' @import shiny
  #' 
  #' 
  if(angle <= 90){
    adj.ang <- (angle*pi)/180
    y.shift <- dista * cos(adj.ang)
    x.shift <- dista * sin(adj.ang)
    point.coords <- c(xcord + x.shift, ycord + y.shift)
  }else{
    if(angle > 90 & angle <= 180){
      adj.ang <-  ((angle - 90)*pi)/180
      x.shift <- dista * cos(adj.ang)
      y.shift <- dista * sin(adj.ang)
      point.coords <- c(xcord + x.shift, ycord - y.shift)
    }else{
      if(angle > 180 & angle <= 270){
        adj.ang <- ((angle - 180)*pi)/180 
        x.shift <- dista * sin(adj.ang)
        y.shift <- dista * cos(adj.ang)
        point.coords <- c(xcord - x.shift, ycord - y.shift)
      }else{
        adj.ang <- ((angle -270)*pi)/180
        x.shift <- dista * cos(adj.ang)
        y.shift <- dista * sin(adj.ang)
        point.coords <- c(xcord - x.shift, ycord + y.shift)
      }
    }
  } 
  return(point.coords)
}
