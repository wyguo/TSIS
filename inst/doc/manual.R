## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()

## ----eval=F--------------------------------------------------------------
#  install.packages(c("shiny", "shinythemes","ggplot2","plotly","zoo","gtools"), dependencies=TRUE)
#  

## ----eval=F--------------------------------------------------------------
#  ##if devtools is not installed, typing
#  #install.packages('devtools')
#  
#  library(devtools)
#  devtools::install_github("wyguo/TSIS")

## ----eval=T--------------------------------------------------------------
library(TSIS)

## ----eval=F--------------------------------------------------------------
#  R.home()

## ----echo=T--------------------------------------------------------------
##26 time points, 3 biological replicates and 3 technical replicates, in total 234 sample points. 
AtRTD2$data.exp[1:10,1:3]
AtRTD2$mapping[1:10,]
colnames(AtRTD2$data.exp)[1:10]

## ----echo=T,eval=F-------------------------------------------------------
#  AtRTD2.example(dir='data')

## ------------------------------------------------------------------------
##use function TSIS::rowmean to take average values
data.exp.mean<-rowmean(t(AtRTD2$data.exp),group = paste0('T',rep(1:26,each=9)))
data.exp.mean<-t(data.exp.mean)
data.exp.mean[1:10,1:4]

##example, to find the intersection points of iso1 and iso2
iso1='AT1G13350_ID2'
iso2='AT1G13350_P1'

##x1 and x2 are the numeric values for two isoforms to search for the intersection points
##x.points and y.points are the x axis and y axis coordinate values for the time course intersection points. 
ts.intersection(x1=as.numeric(data.exp.mean[iso1,]),x2=as.numeric(data.exp.mean[iso2,]))


## ----echo=T--------------------------------------------------------------
##use function TSIS::ts.spline to fit the samples with smooth curve.
##estimate the values at time points 1-26 on the fitted curves
data.exp.splined<-apply(AtRTD2$data.exp[1:10,],1,
                        function(x) ts.spline(x,t.start = 1,t.end = 26,nrep = 9,df = 18))
data.exp.splined<-t(data.exp.splined)
data.exp.splined[,1:4]

##x1 and x2 are the numeric values for two isoforms to search for the intersection points
##x.points and y.points are the x axis and y axis coordinate values for the time course intersection points. 
ts.intersection(x1=as.numeric(data.exp.splined[iso1,]),x2=as.numeric(data.exp.splined[iso2,]))

## ------------------------------------------------------------------------
t.test(c(1,1,2,2,3,4),c(3,4,5,5,6,6),paired = T)$p.value

## ------------------------------------------------------------------------
cor(c(1,2,3,3,5,6,4,5,6,1,2,3),c(4,5,6,1,2,4,1,2,3,4,5,6),method = 'pearson')

## ------------------------------------------------------------------------
##load the data
data.exp<-AtRTD2$data.exp
mapping<-AtRTD2$mapping
dim(data.exp)
dim(mapping)


## ------------------------------------------------------------------------
##Scores
scores.mean2int<-iso.switch(data.exp=data.exp,mapping =mapping,
                     t.start=1,t.end=26,nrep=9,rank=F,
                     min.t.points =2,min.distance=1,spline =F,spline.df = 9,verbose = F)

## ------------------------------------------------------------------------
##Scores
scores.spline2int<-suppressMessages(iso.switch(data.exp=data.exp,mapping =mapping,
                     t.start=1,t.end=26,nrep=9,rank=F,
                     min.t.points =2,min.distance=1,spline =T,spline.df = 9,verbose = F))

## ------------------------------------------------------------------------
##intersection from mean expression
scores.mean2int.filtered<-score.filter(scores = scores.mean2int,prob.cutoff = 0.5,dist.cutoff = 1,
                                       t.points.cutoff = 2,pval.cutoff = 0.01, cor.cutoff = 0.5,
                                       data.exp = NULL,mapping = NULL,sub.isoform.list = NULL,
                                       sub.isoform = F,max.ratio = F,x.value.limit = c(9,17) )

scores.mean2int.filtered[1:5,]

##intersection from spline method
scores.spline2int.filtered<-score.filter(scores = scores.spline2int,prob.cutoff = 0.5,
                                         dist.cutoff = 1,t.points.cutoff = 2,pval.cutoff = 0.01,
                                         cor.cutoff = 0.5,data.exp = NULL,mapping = NULL,
                                         sub.isoform.list = NULL,sub.isoform = F,max.ratio = F,
                                         x.value.limit = c(9,17) )
  

## ------------------------------------------------------------------------
##intersection from mean expression
sub.isoform.list<-AtRTD2$sub.isoforms
sub.isoform.list[1:10]
scores.mean2int.filtered.subset<-score.filter(scores = scores.mean2int,prob.cutoff = 0.5,dist.cutoff = 1,
                                       t.points.cutoff = 2,pval.cutoff = 0.01, cor.cutoff = 0.5,
                                       data.exp = NULL,mapping = NULL,sub.isoform.list = sub.isoform.list,
                                       sub.isoform = T,max.ratio = F,x.value.limit = c(9,17) )

scores.mean2int.filtered.subset[1:5,]


## ----fig.width=8.5,fig.height=4------------------------------------------
plotTSIS(data2plot = data.exp,scores = scores.mean2int.filtered,iso1 = 'AT3G61600_P1',
        iso2 = 'AT3G61600_P2',gene.name = NULL,y.lab = 'Expression',make.plotly = F,
        t.start = 1,t.end = 26,nrep = 9,prob.cutoff = 0.5,x.lower.boundary = 9,
        x.upper.boundary = 17,show.region = T,show.scores = T,
        line.width =0.5,point.size = 3,error.type = 'stderr',show.errorbar = T,errorbar.size = 0.5,
        errorbar.width = 0.2,spline = F,spline.df = NULL,ribbon.plot = F )


## ----fig.width=8.5,fig.height=4------------------------------------------
plotTSIS(data2plot = data.exp,scores = scores.mean2int.filtered,iso1 = 'AT3G61600_P1',
        iso2 = 'AT3G61600_P2',gene.name = NULL,y.lab = 'Expression',make.plotly = F,
        t.start = 1,t.end = 26,nrep = 9,prob.cutoff = 0.5,x.lower.boundary = 9,
        x.upper.boundary = 17,show.region = T,show.scores = T,error.type = 'stderr',
        line.width =0.5,point.size = 3,show.errorbar = T,errorbar.size = 0.5,
        errorbar.width = 0.2,spline = F,spline.df = NULL,ribbon.plot = T )


## ----eval=F--------------------------------------------------------------
#  TSIS.app(data.size.max = 100)

## ----session, echo=FALSE-------------------------------------------------
sessionInfo()

