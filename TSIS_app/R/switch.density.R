#' Density plot
#'
#' @param x a numeric vector to plot th density.
#' @param time.points a numeric vector of the time points of time-series, e.g. 1,2,3,...
#' @param make.plotly logical, to plot \code{\link{plotly}} format figures (TRUE)
#' or general plot (FALSE)?. See details in \code{\link{ggplotly}} in \code{\link{plotly}} R pacakge.
#' @param plot.type the plot types. Options are "density" for density bar plot and "frequency" for frequency bar plot.
#' @param show.line logical, to show density or frequency line or not?
#' @param titile the title of the plot.
#' @param ... additional parameters pass to \code{\link{plotly::ggplotly}}
#'
#' @return density plot in ggplot2 format or \code{\link{plotly}} format if \code{make.plotly=T}
#'
#' @export
#'



switch.density<-function(x,time.points,
                         make.plotly=T,
                         plot.type='density',
                         show.line=T,
                         title="Density of switch points",...){

  ##plot type
  plot.type<-match.arg(plot.type,c('density','frequency'))

  data2plot<-data.frame(x=x+0.001)

  if(plot.type=='density'){
    ##density plot
    g<-ggplot(data2plot,aes(x,fill=1))+geom_histogram(aes(y=..density..),breaks=time.points,closed='left')+theme_bw()+theme(legend.position='none')
    if(show.line)
      g<-g+stat_density(geom='line',size=1,color='red')
  } else {
    ##frequency plot
    g<-ggplot(data2plot,aes(x,fill=1))+geom_histogram(aes(y=..count..),breaks=time.points,closed='left')+theme_bw()+theme(legend.position='none')
    if(show.line)
      g<-g+geom_freqpoly(binwidth=1,color='red',size=1)
  }

  g<-g+coord_cartesian(xlim=c(min(time.points),max(time.points)))+labs(x='Switch time points',title=title)

  if(make.plotly)
    plotly::ggplotly(g,...) else g
}


