#' Fit a time-series with spline curve
#'
#' @param x a vector of time-series expression
#' @param t.start,t.end start and end time points of the time-series data. The time step is assumed to be 1.
#' @param nrep number of replicates.
#' @param df degree of freedom used in \code{\link{ns}} in \code{\link{splines}} package.
#' @param ... additional arguments passed to \code{\link{predict}}.
#'
#' @return predictions results returned from \code{\link{predict}}.
#'
#' @seealso \code{\link{lm}}, \code{\link{ns}}, \code{\link{predict}}
#'
#' @export
#'

ts.spline<-function(x,t.start,t.end,nrep,df=5,...){
  df<-min(df,t.end-t.start-1)
  times<-rep(t.start:t.end,each=nrep)
  x<-data.frame(times=times,value=as.numeric(x))
  fit<-lm(value~splines::ns(times,df=df),data=x)
  predict(fit,data.frame(times=t.start:t.end),...)
}

