#' Density plot
#'
#' @param x a numeric vector to plot th density.
#' @param t.start,t.end start and end time points. Time step is assumed to be 1.
#' @param make.plotly logical, to plot \code{\link{plotly}} format figures (TRUE)
#' or general plot (FALSE)?. See details in \code{\link{ggplotly}} in \code{\link{plotly}} R pacakge.
#' @param plot.type the plot types. Options are "Density_line" for density line plot, "Density_bar" for density bar plot and "Frequency_bar" for frequency bar plot.
#' @param ... additional parameters pass to \code{\link{density}}
#'
#' @return density plot in ggplot2 format or \code{\link{plotly}} format if \code{make.plotly=T}
#'
#' @export
#'



switch.density<-function(x,t.start=1,t.end=26,make.plotly=T,plot.type='density',show.line=T,title="Density of switch points",...){
  ##create density/frequency lines
  if(show.line)
    y.line<-hist(x,breaks = 1000*(t.end-t.start),plot = F)

  ##plot type
  plot.type<-match.arg(plot.type,c('density','frequency'))

  data2plot<-data.frame(x=x)

  if(plot.type=='density'){
    ##density plot
    g<-ggplot(data2plot,aes(x,fill=1))+geom_histogram(aes(y=..density..),binwidth=1)+theme_bw()+theme(legend.position='none')
    if(show.line)
      g<-g+stat_density(geom='line',size=1,color='red')
  } else {
    ##frequency plot
    g<-ggplot(data2plot,aes(x,fill=1))+geom_histogram(aes(y=..count..),binwidth=1)+theme_bw()+theme(legend.position='none')
    if(show.line)
      g<-g+geom_freqpoly(binwidth=1,color='red',size=1)
  }

  g<-g+coord_cartesian(xlim=c(t.start,t.end))+labs(x='Switch time points')

  if(make.plotly)
    plotly::ggplotly(g) else g
}

# switch.density<-function(x,t.start=1,t.end=26,make.plotly=T,plot.type='Density_line',title="Density of switch points",...){
#   if(plot.type=='Density_line'){
#   x<-density(x,...)
#   x.lab <- 'switch time'
#   y.lab<-'Density'
#   if(make.plotly)
#     plotly::plot_ly(data.frame(x=x$x,y=x$y),x=~x,y=~y) %>%
#     add_lines() %>%
#     layout(title = title,xaxis = list(title=x.lab),yaxis=list(title=y.lab))
#   else plot(x,col='blue',main=title,xlab = x.lab,ylab=y.lab)
#   } else {
#     y<-hist(x,breaks = t.start:t.end,plot = F)
#     if(plot.type=='Density_bar')
#       data2plot<-data.frame(x=y$breaks[-1],y=y$density)
#     if(plot.type=='Frequency_bar')
#       data2plot<-data.frame(x=y$breaks[-1],y=y$counts)
#
#     g<-ggplot(data2plot,aes(x=x,y=y,fill=1))+geom_bar(stat = 'identity')+theme_bw()+theme(legend.position='none')+
#       labs(x='switch time',y=strsplit(plot.type,'_')[[1]][1],title=title)
#
#     if(make.plotly)
#       plotly::ggplotly(g)
#     else g
#
#   }
#
#
# }
