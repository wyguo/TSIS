#' Plot time-series switches of a pair of isoforms
#'
#' @param data2plot the expession data frame used to plot the isoforms. The row names of the data frame are the isoforms to plot and the column names are
#' time points and the replicates of time points.The column ordering should be for example t1.replicate1, t1.replicate2, ..., t2.replicate1, t2.replicate2, ...
#' @param scores the output scores of functions \code{\link{iso.switch}} and \code{\link{score.filter}}. Default is NULL. To show score labels on the plot,
#' the scores must be provided.
#' @param iso1,iso2 character string names of the first and second isoforms to plot. If \code{iso1} and \code{iso2} are NULL,
#' the input \code{data2plot} must be a data frame of two rows and the row names are used as isoforms to plot.
#' @param gene.name a character string of gene name to show as the title of the plot. If \code{gene.name=NULL}, the plot title will be "iso1_vs_iso2".
#' @param y.lab the y label of the plot, default is "Expression".
#' @param make.plotly logical, to plot \code{\link{plotly}} format figures (TRUE) or plain ggplot2 format figures(FALSE)?. See details in \code{\link{ggplotly}} in \code{\link{plotly}} R pacakge.
#' @param t.start,t.end start and end time points of the time-series data. The time step is assumed to be 1.
#' @param nrep number of replicates.
#' @param prob.cutoff the cut-off of switch frequencies/probabilities to label the switch points.
#' @param x.lower.boundary,x.upper.boundary the lower and upper boundary of x axis (time points) for the region under investigation. In the analysis, if the isoform pairs have no intersection points
#' in this region, they are filtered.
#' @param show.region logical, to highlight the region under investigation (TRUE) or not (FALSE)?
#' @param show.scores logical, to show score labels on the plot (TRUE) or not (FALSE)? The \code{scores} object must be provided.
#' @param error.type the error type used to show the error bar in the plots. Options are "stderr" for standard error and "sd" for standard deviation. See details in \code{\link{data.error}}.
#' @param show.errorbar logical, to show error bar (TRUE) or not (FALSE) in the error bar plot.
#' @param errorbar.width,errorbar.size the width and size of error bars. See detials in \code{\link{geom_errorbar}} in \code{\link{ggplot2}} R pacakge.
#' @param line.width,point.size line width and point marker size of the plots.
#' @param spline logical, to plot the spline smoothed lines (TRUE) or the lines of mean expression (FALSE).
#' @param spline.df the degree of freedom used for spline.  The value must be the same as in the function \code{\link{iso.switch}}. See spline details in \code{\link{geom_smooth}} in \code{\link{ggplot2}}
#' and \code{\link{ns}} in \code{\link{splines}}.
#' @param ribbon.plot logcial, to make ribbon plot (TRUE) or error bar plot (FALSE). See ribbon plot details in \code{\link{geom_smooth}} in \code{\link{ggplot2}} R pacakge.
#'
#'
#'
#' @examples
#'
#' #generate random datasets
#' set.seed(100)
#' data2plot<-data.frame(matrix(abs(rnorm(60)),nrow=2),row.names = c('iso1','iso2'))
#' #error bar plot
#' plotTSIS(data2plot=data2plot,scores=NULL,iso1=NULL,iso2=NULL,gene.name=NULL,
#'         y.lab='Expression',make.plotly=F,t.start=1,t.end=10,nrep=3,prob.cutoff=0.5,x.lower.boundary=3,
#'         x.upper.boundary=8,show.region=T,error.type='stderr',ribbon.plot = F)
#'
#' #ribbon plot
#' plotTSIS(data2plot=data2plot,scores=NULL,iso1=NULL,iso2=NULL,gene.name=NULL,
#'          y.lab='Expression',make.plotly=F,t.start=1,t.end=10,nrep=3,prob.cutoff=0.5,x.lower.boundary=3,
#'          x.upper.boundary=8,show.region=T,error.type='stderr',ribbon.plot = T)
#'
#' @return a ggplot2 or \code{\link{ggplotly}} (if \code{plotly=TRUE}) plot.
#' @export
#'
#' @seealso \code{\link{ggplotly}}, \code{\link{iso.switch}}, \code{\link{score.filter}}, \code{\link{data.error}}, \code{\link{geom_smooth}}, \code{\link{ns}}
#'
plotTSIS<-function(data2plot,scores=NULL,iso1=NULL,iso2=NULL,gene.name=NULL,y.lab='Expression',make.plotly=F,
                  t.start=1,t.end=26,nrep=9,prob.cutoff=0.5,x.lower.boundary=9,x.upper.boundary=17,show.region=T,show.scores=T,
                  error.type='stderr',show.errorbar=T,errorbar.width=0.2,errorbar.size=0.5,line.width=1,point.size=3,spline=F,spline.df=NULL,ribbon.plot=F){

  require(ggplot2)

  if(is.null(spline.df))
    spline.df<-floor((t.end-t.start)*2/3)

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  y.upper <- function(x){mean(x) + data.error(x)}
  y.lower<-function(x){mean(x) - data.error(x)}


  if(is.null(iso1) | is.null(iso2)){
    iso1<-rownames(data2plot)[1]
    iso2<-rownames(data2plot)[2]
  }

  data2plot<-data2plot[c(iso1,iso2),]

  if(is.null(gene.name))
    gene.name<-paste0(iso1,'_vs_',iso2)
  ##arange the input data to plot data
  colnames(data2plot)<-paste0(rep(t.start:t.end,each=nrep),'_',rep(1:nrep,(t.end-t.start+1)))
  data2plot<-reshape2::melt(as.matrix(data2plot))
  idx<-do.call(rbind,strsplit(as.vector(data2plot[,2]),split = '_'))
  data2plot<-data.frame(isoforms=data2plot[,1],times=as.numeric(idx[,1]),reps=as.numeric(idx[,2]),value=data2plot[,3])

  ##ribbon plot
  if(ribbon.plot){
    if(spline){
      data2plot<-data.frame(data2plot[,1:2],mean=data2plot$value,error=0)
      data2plot$times<-factor(data2plot$times,levels = unique(data2plot$times))
      g<-ggplot(data2plot, aes(x=times, y=mean,group=isoforms,color=isoforms,fill=isoforms,shape=isoforms)) +
        geom_point() +theme_bw()+
        geom_smooth(method = "lm",formula = y ~ splines::ns(x, spline.df),size=line.width)
    } else {
      data2plot<-data.frame(data2plot[,1:2],mean=data2plot$value,error=0)
      data2plot$times<-factor(data2plot$times,levels = unique(data2plot$times))
      g=ggplot(data2plot, aes(x=times, y=mean,group=isoforms,color=isoforms,fill=isoforms,shape=isoforms)) +
        geom_point() +theme_bw()+
        stat_summary(fun.y=mean,geom='line',size=line.width)+
        stat_summary(fun.y=mean,geom = 'ribbon',fun.ymax = y.upper,fun.ymin = y.lower,alpha=0.3,colour=NA)
    }
  ##error bar plots
  } else{
    if(spline){
      values<-by(data2plot$value,INDICES = data2plot$isoforms,simplify = T,
                 FUN = function(x) ts.spline(x,t.start = t.start,t.end = t.end,nrep = nrep,df=spline.df,se.fit=T))
      x1=data.frame(mean.Group.1=paste0(iso1,'_at_',t.start:t.end),mean.x=values[[iso1]]$fit,error=values[[iso1]]$se.fit)
      x2=data.frame(mean.Group.1=paste0(iso2,'_at_',t.start:t.end),mean.x=values[[iso2]]$fit,error=values[[iso2]]$se.fit)
      data2plot<-rbind(x1,x2)
    } else {
      data2plot<-with(data2plot,{
        data.frame(
          mean=aggregate(value,by = list(interaction(isoforms,times,sep = '_at_')),mean),
          error=aggregate(value,by = list(interaction(isoforms,times,sep = '_at_')),function(x) data.error(x,error.type = error.type))[,-1]
        )
      })

    }


    data2plot<-data.frame(do.call(rbind,strsplit(as.vector(data2plot[,1]),split = '_at_')),data2plot[,-1])
    colnames(data2plot)<-c('isoforms','times','mean','error')
    data2plot$times<-factor(data2plot$times,levels = unique(data2plot$times))

    g<-ggplot(data2plot,aes(x=times,y=mean,group=isoforms,shape=isoforms,color=isoforms))+theme_bw()+
      geom_line(size=line.width)+geom_point(size=point.size,aes(fill=isoforms))
    if(show.errorbar)
      g<-g+geom_errorbar(data=data2plot,aes(ymin=mean-error,ymax=mean+error),width=errorbar.width,color='black',size=errorbar.size)
  }



  if(is.null(scores))
    idx=NULL else idx<-which(scores$iso1 %in% c(iso1,iso2) & scores$iso2 %in% c(iso1,iso2) & scores$prob>=prob.cutoff)

  if(length(idx)==0){
    g<-g } else {
      sub.scores<-scores[idx,]
      data2points<-data.frame(isoforms='switch_points',times=sub.scores$x.value,mean=sub.scores$y.value,error=0)

      g<-g+geom_point(data=data2points,aes(x=times,y=mean,color=isoforms,fill=isoforms,shape=isoforms),size=point.size)+
        scale_color_manual(values=c(gg_color_hue(2),rep('black',length(idx))),breaks=c(iso1,iso2,rep('switch_points',length(idx))))+
        scale_fill_manual(values=c(gg_color_hue(2),rep('black',length(idx))),breaks=c(iso1,iso2,rep('switch_points',length(idx))))+
        scale_shape_manual(values=c(15,17,rep(16,length(idx))))
      if(show.scores)
        g<-g+annotate('text',x = 1.1*sub.scores$x.value, y =sub.scores$y.value+max(data2plot$mean)/20,
                      label = paste0('prob=',round(sub.scores$prob,2),'; dist=',round(sub.scores$dist,2),'; cor=',round(sub.scores$cor,2)))
    }

  if(show.region){
    g<-g+annotate("rect", xmin = x.lower.boundary, xmax = x.upper.boundary, ymin = -Inf, ymax = Inf,alpha = 0.2,color='skyblue',fill='skyblue')+
      annotate("text",x=mean(c(x.upper.boundary,x.lower.boundary)),
               y=max(data2plot$mean+data2plot$error),label='Region for investigation',color='blue',size=5)+
      geom_vline(xintercept = c(x.lower.boundary,x.upper.boundary),linetype=3,lwd=1,color='skyblue',alpha=1)
  }
  g<-g+labs(title=paste("Isoforms of gene:",gene.name),x='Time points',y=y.lab)

  if(make.plotly)
    plotly::ggplotly(g) else g
}


