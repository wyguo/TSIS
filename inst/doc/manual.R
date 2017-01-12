## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()

## ----eval=F--------------------------------------------------------------
#  install.packages(c("shiny", "shinythemes","ggplot2","plotly","zoo","gtools","devtools"), dependencies=TRUE)
#  

## ----eval=F--------------------------------------------------------------
#  library(devtools)
#  devtools::install_github("wyguo/TSIS")

## ----eval=T--------------------------------------------------------------
library(TSIS)

## ----eval=F--------------------------------------------------------------
#  R.home()

## ----eval=F--------------------------------------------------------------
#  ##26 time points, 3 biological replicates and 3 technical replicates, in total 234 sample points.
#  library(TSIS)
#  AtRTD2$data.exp[1:10,1:3]
#  AtRTD2$mapping[1:10,]
#  AtRTD2$sub.isoforms[1:10]

## ----echo=T,eval=F-------------------------------------------------------
#  AtRTD2.example()

## ----echo=T,eval=T-------------------------------------------------------
##load the data
data.exp<-AtRTD2$data.exp
mapping<-AtRTD2$mapping
dim(data.exp);dim(mapping)

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



## ----eval=F--------------------------------------------------------------
#  
#  ##intersection from spline method
#  scores.spline2int.filtered<-score.filter(scores = scores.spline2int,prob.cutoff = 0.5,
#                                           dist.cutoff = 1,t.points.cutoff = 2,pval.cutoff = 0.01,
#                                           cor.cutoff = 0.5,data.exp = NULL,mapping = NULL,
#                                           sub.isoform.list = NULL,sub.isoform = F,max.ratio = F,
#                                           x.value.limit = c(9,17) )
#  

## ----eval=F--------------------------------------------------------------
#  ##intersection from mean expression
#  sub.isoform.list<-AtRTD2$sub.isoforms
#  sub.isoform.list[1:10]
#  scores.mean2int.filtered.subset<-score.filter(scores = scores.mean2int,prob.cutoff = 0.5,dist.cutoff = 1,
#                                         t.points.cutoff = 2,pval.cutoff = 0.01, cor.cutoff = 0.5,
#                                         data.exp = NULL,mapping = NULL,sub.isoform.list = sub.isoform.list,
#                                         sub.isoform = T,max.ratio = F,x.value.limit = c(9,17) )
#  

## ----eval=F--------------------------------------------------------------
#  scores.mean2int.filtered.maxratio<-score.filter(scores = scores.mean2int,prob.cutoff = 0.5,dist.cutoff = 1,
#                                         t.points.cutoff = 2,pval.cutoff = 0.01, cor.cutoff = 0,
#                                         data.exp = data.exp,mapping = mapping,sub.isoform.list = NULL,
#                                         sub.isoform = F,max.ratio = T,x.value.limit = c(9,17) )

## ----eval=F,fig.width=8.5,fig.height=4-----------------------------------
#  plotTSIS(data2plot = data.exp,scores = scores.mean2int.filtered,iso1 = 'AT3G61600_P1',
#          iso2 = 'AT3G61600_P2',gene.name = NULL,y.lab = 'Expression',make.plotly = F,
#          t.start = 1,t.end = 26,nrep = 9,prob.cutoff = 0.5,x.lower.boundary = 9,
#          x.upper.boundary = 17,show.region = T,show.scores = T,
#          line.width =0.5,point.size = 3,error.type = 'stderr',show.errorbar = T,errorbar.size = 0.5,
#          errorbar.width = 0.2,spline = F,spline.df = NULL,ribbon.plot = F )
#  

## ----eval=F,fig.width=8.5,fig.height=4-----------------------------------
#  plotTSIS(data2plot = data.exp,scores = scores.mean2int.filtered,iso1 = 'AT3G61600_P1',
#          iso2 = 'AT3G61600_P2',gene.name = NULL,y.lab = 'Expression',make.plotly = F,
#          t.start = 1,t.end = 26,nrep = 9,prob.cutoff = 0.5,x.lower.boundary = 9,
#          x.upper.boundary = 17,show.region = T,show.scores = T,error.type = 'stderr',
#          line.width =0.5,point.size = 3,show.errorbar = T,errorbar.size = 0.5,
#          errorbar.width = 0.2,spline = F,spline.df = NULL,ribbon.plot = T )
#  

## ----session, echo=FALSE-------------------------------------------------
sessionInfo()

