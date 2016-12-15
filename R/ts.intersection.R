#' Searching for intersection points of two time-series
#'
#' @param x1 a vector of values for first time-series
#' @param x2 a vector of values for second time-series
#'
#' @return a data frame of intersection points. The first and second columns are x and y coordinates values of intersections, respectively.
#'
#' @export
#'
ts.intersection<-function(x1,x2){
  # Points always intersect when above=TRUE, then FALSE or reverse
  # Find points where x1 is above x2.
  above<-x1>x2
  intersect.points<-which(diff(above)!=0)
  if(length(intersect.points)==0)
    return(NULL)
  # Find the slopes for each line segment.
  x1.slopes<-x1[intersect.points+1]-x1[intersect.points]
  x2.slopes<-x2[intersect.points+1]-x2[intersect.points]
  # Find the intersection for each segment.
  x.points<-intersect.points + ((x2[intersect.points] - x1[intersect.points]) / (x1.slopes-x2.slopes))
  y.points<-x1[intersect.points] + (x1.slopes*(x.points-intersect.points))
  ##intersection points
  intersection.points=data.frame(x.points=x.points,y.points=y.points,row.names = NULL)
  intersection.points<-intersection.points[!duplicated(intersection.points),]
  return(intersection.points)
}
