## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()

## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache = TRUE,eval=F,echo = T)
figure.type='word'

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(figure.type=='html')
  cat('<img src="fig/manual figures_001.png" width="900px" height="200px" />')
if(figure.type=='word')
  cat('![](fig/manual figures_001.png)')
if(figure.type=='manual')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_001.png)') 

## ------------------------------------------------------------------------
#  R.home()

## ------------------------------------------------------------------------
#  library(devtools)
#  devtools::install_github("wyguo/TSIS")

## ------------------------------------------------------------------------
#  library(TSIS)

## ------------------------------------------------------------------------
#  TSIS.app()

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(figure.type=='html')
  cat('<img src="fig/manual figures_002.png" width="600px" height="50px" />')
if(figure.type=='word')
  cat('![](fig/manual figures_002.png)')
if(figure.type=='manual')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_002.png)') 

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(figure.type=='html')
  cat('<img src="fig/manual figures_003.png" width="800px" height="400px" />')
if(figure.type=='word')
  cat('![](fig/manual figures_003.png)')
if(figure.type=='manual')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_003.png)') 

## ------------------------------------------------------------------------
#  ##26 time points, 3 biological replicates and 3 technical replicates, in total 234 sample points.
#  library(TSIS)
#  AtRTD2$data.exp[1:10,1:3]
#  AtRTD2$mapping[1:10,]
#  AtRTD2$sub.isoforms[1:10]

## ------------------------------------------------------------------------
#  AtRTD2.example()

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(figure.type=='html')
  cat('<img src="fig/manual figures_004.png" width="600px" height="300px" />')
if(figure.type=='word')
  cat('![](fig/manual figures_004.png)')
if(figure.type=='manual')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_004.png)') 

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(figure.type=='html')
  cat('<img src="fig/manual figures_005.png" width="800px" height="400px" />')
if(figure.type=='word')
  cat('![](fig/manual figures_005.png)')
if(figure.type=='manual')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_005.png)') 

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(figure.type=='html')
  cat('<img src="fig/manual figures_006.png" width="800px" height="400px" />')
if(figure.type=='word')
  cat('![](fig/manual figures_006.png)')
if(figure.type=='manual')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_006.png)') 

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(figure.type=='html')
  cat('<img src="fig/manual figures_007.png" width="800px" height="400px" />')
if(figure.type=='word')
  cat('![](fig/manual figures_007.png)')
if(figure.type=='manual')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_007.png)') 

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(figure.type=='html')
  cat('<img src="fig/manual figures_008.png" width="800px" height="400px" />')
if(figure.type=='word')
  cat('![](fig/manual figures_008.png)')
if(figure.type=='manual')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_008.png)') 

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(figure.type=='html')
  cat('<img src="fig/manual figures_009.png" width="800px" height="400px" />')
if(figure.type=='word')
  cat('![](fig/manual figures_009.png)')
if(figure.type=='manual')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_009.png)') 

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(figure.type=='html')
  cat('<img src="fig/manual figures_010.png" width="800px" height="400px" />')
if(figure.type=='word')
  cat('![](fig/manual figures_010.png)')
if(figure.type=='manual')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_010.png)') 

## ----eval=F,echo=T-------------------------------------------------------
#  ##load the data
#  library(TSIS)
#  data.exp<-AtRTD2$data.exp
#  mapping<-AtRTD2$mapping
#  dim(data.exp);dim(mapping)

## ----echo=T,eval=F-------------------------------------------------------
#  ##Scores
#  scores.mean2int<-iso.switch(data.exp=data.exp,mapping =mapping,
#                       t.start=1,t.end=26,nrep=9,rank=F,
#                       min.t.points =2,min.distance=1,spline =F,spline.df = 9,verbose = F)

## ------------------------------------------------------------------------
#  ##Scores, set spline=T and define spline degree of freedom to spline.df=9
#  scores.spline2int<-iso.switch(data.exp=data.exp,mapping =mapping,
#                       t.start=1,t.end=26,nrep=9,rank=F,
#                       min.t.points =2,min.distance=1,spline =T,spline.df = 9,verbose = F)
#  

## ----eval=F--------------------------------------------------------------
#  ##intersection from mean expression
#  scores.mean2int.filtered<-score.filter(scores = scores.mean2int,prob.cutoff = 0.5,dist.cutoff = 1,
#                                         t.points.cutoff = 2,pval.cutoff = 0.01, cor.cutoff = 0.5,
#                                         data.exp = NULL,mapping = NULL,sub.isoform.list = NULL,
#                                         sub.isoform = F,max.ratio = F,x.value.limit = c(9,17) )
#  
#  scores.mean2int.filtered[1:5,]

## ------------------------------------------------------------------------
#  ##intersection from spline method
#  scores.spline2int.filtered<-score.filter(scores = scores.spline2int,prob.cutoff = 0.5,
#                                           dist.cutoff = 1,t.points.cutoff = 2,pval.cutoff = 0.01,
#                                           cor.cutoff = 0.5,data.exp = NULL,mapping = NULL,
#                                           sub.isoform.list = NULL,sub.isoform = F,max.ratio = F,
#                                           x.value.limit = c(9,17) )

## ------------------------------------------------------------------------
#  ##intersection from mean expression
#  ##input a list of isoform names for investigation.
#  sub.isoform.list<-AtRTD2$sub.isoforms
#  sub.isoform.list[1:10]
#  ##assign the isoform name list to sub.isoform.list and set sub.isoform=TRUE
#  scores.mean2int.filtered.subset<-score.filter(scores = scores.mean2int,prob.cutoff = 0.5,dist.cutoff = 1,
#                                         t.points.cutoff = 2,pval.cutoff = 0.01, cor.cutoff = 0.5,
#                                         data.exp = NULL,mapping = NULL,sub.isoform.list = sub.isoform.list,
#                                         sub.isoform = T,max.ratio = F,x.value.limit = c(9,17) )

## ------------------------------------------------------------------------
#  scores.mean2int.filtered.maxratio<-score.filter(scores = scores.mean2int,prob.cutoff = 0.5,dist.cutoff = 1,
#                                         t.points.cutoff = 2,pval.cutoff = 0.01, cor.cutoff = 0,
#                                         data.exp = data.exp,mapping = mapping,sub.isoform.list = NULL,
#                                         sub.isoform = F,max.ratio = T,x.value.limit = c(9,17) )

## ------------------------------------------------------------------------
#  plotTSIS(data2plot = data.exp,scores = scores.mean2int.filtered,iso1 = 'AT5G60930_P2',
#          iso2 = ' AT5G60930_P3',gene.name = NULL,y.lab = 'Expression',make.plotly = F,
#          t.start = 1,t.end = 26,nrep = 9,prob.cutoff = 0.5,x.lower.boundary = 9,
#          x.upper.boundary = 17,show.region = T,show.scores = T,
#          line.width =0.5,point.size = 3,error.type = 'stderr',show.errorbar = T,errorbar.size = 0.5,
#          errorbar.width = 0.2,spline = F,spline.df = NULL,ribbon.plot = F )

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(figure.type=='html')
  cat('<img src="fig/AT5G60930_P2_1.png" width="800px" height="400px" />')
if(figure.type=='word')
  cat('![](fig/AT5G60930_P2_1.png)')
if(figure.type=='manual')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/AT5G60930_P2_1.png)') 

## ------------------------------------------------------------------------
#  plotTSIS(data2plot = data.exp,scores = scores.mean2int.filtered,iso1 = 'AT5G60930_P2',
#          iso2 = 'AT5G60930_P3',gene.name = NULL,y.lab = 'Expression',make.plotly = F,
#          t.start = 1,t.end = 26,nrep = 9,prob.cutoff = 0.5,x.lower.boundary = 9,
#          x.upper.boundary = 17,show.region = T,show.scores = T,error.type = 'stderr',
#          line.width =0.5,point.size = 3,show.errorbar = T,errorbar.size = 0.5,
#          errorbar.width = 0.2,spline = F,spline.df = NULL,ribbon.plot = T )
#  

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(figure.type=='html')
  cat('<img src="fig/AT5G60930_P2_2.png" width="800px" height="400px" />')
if(figure.type=='word')
  cat('![](fig/AT5G60930_P2_2.png)')
if(figure.type=='manual')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/AT5G60930_P2_2.png)') 

## ----session, echo=FALSE,eval=T------------------------------------------
sessionInfo()

