#' Calculate column mean of a matrix or data frame based on a grouping variable
#'
#' Compute column means across rows of a numeric matrix-like object for each level
#' of a grouping variable.
#' @param x a matrix or data frame.
#' @param group a vector of factor giving grouping, with one element per row of x.
#' @param reorder if TRUE, then the result will be in order of sort(unique(group)).
#' @param na.rm logical (TRUE or FALSE). Should NA (including NaN) values be replaced by value 0?
#' @return \code{rowmean} returns a matrix or data frame containing the means. There
#' will be one row per unique value of group.
#'
#' @examples
#' x <- matrix(runif(50), ncol = 5)
#' group <- sample(1:4, 10, TRUE)
#' xmean <- rowmean(x, group)
#'
#' @export
#' @seealso \code{\link{rowsum}}, \code{\link{rowratio}}
rowmean<-function(x,group,reorder=T,na.rm=T){
  # counts<-as.vector(t(rowsum(data.frame(x=rep(1,nrow(x))),group = group)))
  counts<-as.numeric(table(group))
  sums<-rowsum(x,group = group)
  means<-(diag(1/counts))%*% as.matrix(sums)
  rownames(means)<-unique(group)
  if(na.rm)
    means[is.na(means)]<-0
  if(reorder)
    means<-means[gtools::mixedsort(rownames(means)),]
  return(means)
}



