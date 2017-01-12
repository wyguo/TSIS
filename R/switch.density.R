#' Density plot
#'
#' @param x a numeric vector to plot th density.
#' @param make.plotly logical, to plot \code{\link{plotly}} format figures (TRUE)
#' or general plot (FALSE)?. See details in \code{\link{ggplotly}} in \code{\link{plotly}} R pacakge.
#' @param ... additional parameters pass to \code{\link{density}}
#'
#' @return density plot in general format or \code{\link{plotly}} format if \code{make.plotly=T}
#'
#' @export
#'
#' @examples
#' ##random values
#' set.seed(100)
#' x<-rnorm(100,mean=0,sd=1)
#' ##geneal format
#' switch.density(x,make.plotly =F)
#' ##plotly format
#' switch.density(x,make.plotly =T)
#'
#'
switch.density<-function(x,make.plotly=T,title="Density of switch points",...){
  x<-density(x,...)
  x.lab <- 'switch time'
  y.lab<-'Density'
  if(make.plotly)
    plot_ly(data.frame(x=x$x,y=x$y),x=~x,y=~y) %>%
    add_lines() %>%
    layout(title = title,xaxis = list(title=x.lab),yaxis=list(title=y.lab))
  else plot(x,col='blue',main=title,xlab = x.lab,ylab=y.lab)

}
