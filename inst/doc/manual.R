## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()

## ----setup, include=T,eval=T,echo=F--------------------------------------
knitr::opts_chunk$set(cache = TRUE,eval=F,echo = T)
file.type='html'

## ----results='asis',echo=F,eval=TRUE-------------------------------------
if(file.type=='html')
  cat('<img src="fig/figures_010.png" width="900px"  />')

if(file.type=='word')
  cat('![](fig/figures_010.png)')

if(file.type=='readme')
  cat('![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_010.png)')

