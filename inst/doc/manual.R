## ----setup, include=F,eval=T,echo=F--------------------------------------
knitr::opts_chunk$set(cache = TRUE,eval=F,echo = T)
file.type='html'

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(file.type=='html')
  cat('<img src="fig/figures_001.png" width="900px />')

if(file.type=='word')
  cat('![](fig/figures_001.png)')

if(file.type=='readme')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_001.png)')

## ------------------------------------------------------------------------
#  R.home()

## ------------------------------------------------------------------------
#  library(devtools)
#  devtools::install_github("wyguo/TSIS")

## ------------------------------------------------------------------------
#  library(TSIS)

## ------------------------------------------------------------------------
#  TSIS.data.example()

## ------------------------------------------------------------------------
#  TSIS.app(data.size.max = 100)

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(file.type=='html')
  cat('<img src="fig/figures_002.png" width="900px"  />')

if(file.type=='word')
  cat('![](fig/figures_002.png)')

if(file.type=='readme')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_002.png)')

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(file.type=='html')
  cat('<img src="fig/figures_003.png" width="900px"  />')

if(file.type=='word')
  cat('![](fig/figures_003.png)')

if(file.type=='readme')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_003.png)')

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(file.type=='html')
  cat('<img src="fig/figures_004.png" width="900px"  />')

if(file.type=='word')
  cat('![](fig/figures_004.png)')

if(file.type=='readme')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_004.png)')

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(file.type=='html')
  cat('<img src="fig/figures_005.png" width="900px"  />')

if(file.type=='word')
  cat('![](fig/figures_005.png)')

if(file.type=='readme')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_005.png)')

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(file.type=='html')
  cat('<img src="fig/figures_006.png" width="900px"  />')

if(file.type=='word')
  cat('![](fig/figures_006.png)')

if(file.type=='readme')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_006.png)')

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(file.type=='html')
  cat('<img src="fig/figures_007.png" width="900px"  />')

if(file.type=='word')
  cat('![](fig/figures_007.png)')

if(file.type=='readme')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_007.png)')

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(file.type=='html')
  cat('<img src="fig/figures_008.png" width="900px"  />')

if(file.type=='word')
  cat('![](fig/figures_008.png)')

if(file.type=='readme')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_008.png)')

## ----eval=F,echo=T-------------------------------------------------------
#  ##load the data
#  library(TSIS)
#  data.exp<-TSIS.data$data.exp
#  mapping<-TSIS.data$mapping
#  dim(data.exp);dim(mapping)

## ----echo=T,eval=F-------------------------------------------------------
#  ##Scores
#  scores.mean2int<-iso.switch(data.exp=data.exp,mapping =mapping,
#                       times=rep(1:26,each=9),rank=F,
#                       min.t.points =2,min.difference=1,spline =F,
#                       spline.df = 9,verbose = F)

## ------------------------------------------------------------------------
#  ##Scores, set spline=T and define spline degree of freedom to spline.df=9
#  scores.spline2int<-iso.switch(data.exp=data.exp,mapping =mapping,
#                       times=rep(1:26,each=9),rank=F,
#                       min.t.points =2,min.difference=1,spline =T,
#                       spline.df = 10,verbose = F)

## ----eval=F--------------------------------------------------------------
#  ##intersection from mean expression
#  scores.mean2int.filtered<-score.filter(
#    scores = scores.mean2int,prob.cutoff = 0.5,diff.cutoff = 1,
#    t.points.cutoff = 2,pval.cutoff = 0.001, cor.cutoff = 0,
#    data.exp = NULL,mapping = NULL,sub.isoform.list = NULL,
#    sub.isoform = F,max.ratio = F,x.value.limit = c(1,26)
#  )
#  scores.mean2int.filtered[1:5,]

## ------------------------------------------------------------------------
#  ##intersection from spline method
#  scores.spline2int.filtered<-score.filter(
#    scores = scores.spline2int,prob.cutoff = 0.5,
#    diff.cutoff = 1,t.points.cutoff = 2,pval.cutoff = 0.01,
#    cor.cutoff = 0.5,data.exp = NULL,mapping = NULL,
#    sub.isoform.list = NULL,sub.isoform = F,max.ratio = F,
#    x.value.limit = c(9,17)
#  ))

## ------------------------------------------------------------------------
#  ##intersection from mean expression
#  ##input a list of isoform names for investigation.
#  sub.isoform.list<-TSIS.data$sub.isoforms
#  sub.isoform.list[1:10]
#  ##assign the isoform name list to sub.isoform.list and set sub.isoform=TRUE
#  scores.mean2int.filtered.subset<-score.filter(
#    scores = scores.mean2int,prob.cutoff = 0.5,diff.cutoff = 1,
#    t.points.cutoff = 2,pval.cutoff = 0.01, cor.cutoff = 0.5,
#    data.exp = NULL,mapping = NULL,sub.isoform.list = sub.isoform.list,
#    sub.isoform = T,max.ratio = F,x.value.limit = c(9,17)
#  )

## ------------------------------------------------------------------------
#  scores.mean2int.filtered.maxratio<-score.filter(
#    scores = scores.mean2int,prob.cutoff = 0.5,diff.cutoff = 1,
#    t.points.cutoff = 2,pval.cutoff = 0.01, cor.cutoff = 0,
#    data.exp = data.exp,mapping = mapping,sub.isoform.list = NULL,
#    sub.isoform = F,max.ratio = T,x.value.limit = c(9,17)
#  )

## ------------------------------------------------------------------------
#  library(gridExtra)
#  g1<-switch.density(scores.mean2int.filtered$x.value,make.plotly = F,
#                     show.line = F,plot.type = 'frequency',
#                     title = 'Frequency of switch time' ,time.points = 1:26)
#  g2<-switch.density(scores.mean2int.filtered$x.value,make.plotly = F,
#                     show.line = T,plot.type = 'density',
#                     title = 'Density of switch time' ,time.points = 1:26)
#  gridExtra::grid.arrange(g1,g2,ncol=2)

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(file.type=='html')
  cat('<img src="fig/figures_009.png" width="900px"  />')

if(file.type=='word')
  cat('![](fig/figures_009.png)')

if(file.type=='readme')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_009.png)')

## ------------------------------------------------------------------------
#  plotTSIS(data2plot = data.exp,scores = scores.mean2int.filtered,
#           iso1 = 'G30.iso2',iso2 = 'G30.iso3',gene.name = NULL,
#           y.lab = 'Expression',make.plotly = F,
#           times=rep(1:26,each=9),prob.cutoff = 0.5,
#           x.lower.boundary = 9,x.upper.boundary = 17,
#           show.region = T,show.scores = T,
#           line.width =0.5,point.size = 3,
#           error.type = 'stderr',show.errorbar = T,errorbar.size = 0.5,
#           errorbar.width = 0.2,spline = F,spline.df = NULL,ribbon.plot = F)

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(file.type=='html')
  cat('<img src="fig/figures_010.png" width="900px"  />')

if(file.type=='word')
  cat('![](fig/figures_010.png)')

if(file.type=='readme')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_010.png)')

## ------------------------------------------------------------------------
#  plotTSIS(data2plot = data.exp,scores = scores.mean2int.filtered,
#           iso1 = 'G30.iso2',iso2 = 'G30.iso3',gene.name = NULL,
#           y.lab = 'Expression',make.plotly = F,
#           times=rep(1:26,each=9),prob.cutoff = 0.5,
#           x.lower.boundary = 9,x.upper.boundary = 17,
#           show.region = T,show.scores = T,
#           line.width =0.5,point.size = 3,
#           error.type = 'stderr',show.errorbar = T,errorbar.size = 0.5,
#           errorbar.width = 0.2,spline = F,spline.df = NULL,ribbon.plot = T)

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(file.type=='html')
  cat('<img src="fig/figures_011.png" width="900px"  />')

if(file.type=='word')
  cat('![](fig/figures_011.png)')

if(file.type=='readme')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_010.png)')

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(file.type!='readme')
  cat('#Session Information')

## ----session, echo=FALSE,eval=T------------------------------------------
if(file.type!='readme')
  sessionInfo()

