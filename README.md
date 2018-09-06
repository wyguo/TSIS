---
title: "TSIS: an R package to infer time-series isoform switch of alternative splicing"
subtitle: 'User manual'
author: "Wenbin Guo"
date: "2017-05-24"
output: 
  BiocStyle::html_document
vignette: >
  %\VignetteIndexEntry{Introduction to TSIS}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



# Description
[TSIS](https://github.com/wyguo/TSIS) is an R package for detecting transcript isoform switches in time-series data. Transcript isoform switches occur when a pair of alternatively spliced isoforms reverse the order of their relative expression levels as shown in <a href="#fig1">Figure 1</a>. TSIS characterizes the transcript switch by 1) defining the isoform switch time-points for any pair of transcript isoforms within a gene, 2) describing the switch using five different features or metrics, 3) filtering the results with user’s specifications and 4) visualizing the results using different plots for the user to examine further details of the switches. All the functions are available in the forms of a graphic interface implemented by [Shiny App](https://shiny.rstudio.com/)  (a web application framework for R) ([Chang, et al., 2016](https://shiny.rstudio.com/)), in which users can implement the analysis easily. The tool can also be run using command lines without graphic interface. This tutorial will cover both in the following sections. This manual can be downloaded from: https://github.com/wyguo/TSIS/blob/master/vignettes/tutorial-shiny.zip.

If you use TSIS in your work, please cite: 

Wenbin Guo, Cristiane P. G. Calixto, John W.S. Brown, Runxuan Zhang, "TSIS: an R package to infer alternative splicing isoform switches for time-series data", Bioinformatics, https://doi.org/10.1093/bioinformatics/btx411, 2017. 


<h2 id="fig1"> </h2>

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_001.png)

**Figure 1:**  Isoform switch analysis methods. Expression data with 3 replicates for each condition/time-point is simulated for isoforms $iso_i$  and $iso_j$ (blue and red circles). The points in the plots represent the samples and the black lines connect the average of samples. (A) is the iso-kTSP algorithm for comparisons of two conditions $c_1$  and $c_2$. The Time-Series Isoform Switch (TSIS) tool is designed for detection and characterization of isoform switches for time series data shown in (B). The time-series with 6 time-points is divided into 4 intervals by the intersection points of average expression.

## Determine switch points
Given that a pair of isoforms $iso_i$  and $iso_j$  may have a number of switches in a time-series, we have offered two approaches to search for the switch time-points in TSIS:

-  The first approach takes the average values of the replicates for each time-point for each transcript isoform. Then it searches for the cross points of the average value of two isoforms across the time-points (seen in <a href="#fig1">Figure 1(B)</a>).
- The second approach uses natural spline curves to fit the time-series data for each transcript isoform and find cross points of the fitted curves for each pair of isoforms.

In most cases, these two methods produce very similar results. However, average values of expression may lose precision by not having information of previous and following time-points. The spline method fits the time-series of expression with control points (depending on spline degree of freedom provided) and weights of several neighbours to obtain designed precision (Hastie and Chambers, 1992). The spline method is useful to find global trends in the time-series data when the data is very noisy. However, it may lack details of isoform switch in the local region. It is recommended that users use both average and spline methods to search for the switch points and examine manually when inconsistent results were produced by the above two methods.

## Define switch scoring features

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/Metric.png)
 
# Installation and loading

## Check before installation 
Due to an issue with [devtools](https://cran.r-project.org/web/packages/devtools/index.html), if R software is installed in a directory whose name has a space character in it, e.g. in "C:\\Program Files", users may get error message "'C:\\Program' is not recognized as an internal or external command". This issue has to be solved by making sure that R is installed in a directory whose name has no space characters. Users can check the R installation location by typing

```r
R.home()
```

## Install dependency packages
install.packages(c("shiny", "shinythemes","ggplot2","plotly","zoo","gtools","devtools"), dependencies=TRUE)

## Install TSIS package
Install [TSIS](https://github.com/wyguo/TSIS)  package from Github using [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package.

```r
library(devtools)
devtools::install_github("wyguo/TSIS")
```

## Loading
Once installed, TSIS package can be loaded as normal.

```r
library(TSIS)
```
## Example datasets
The [TSIS](https://github.com/wyguo/TSIS) package provides the example dataset "TSIS.data" with 30 genes and 145 isoforms, analyzed in 26 time-points, each with 9 replicates. The isoform expression is in TPM (transcript per million) format. Other types of transcript quantifications, such as read counts, Percentage Splicing In (PSI) can also be used in [TSIS](https://github.com/wyguo/TSIS).

The data loaded into the [TSIS](https://github.com/wyguo/TSIS) App must be in *.csv format. Users can download the example datasets from https://github.com/wyguo/TSIS/tree/master/data or by typing the following codes in R console:

```r
TSIS.data.example()
```
The data will be saved in a folder "example data" in the working directory. <a href="#fig3">Figure 3</a> shows the examples of input data in csv format. 

# TSIS Shiny App 
To make the implementation more user friendly, TSIS analysis is integrated into a [Shiny App](https://shiny.rstudio.com/) ([Chang, et al., 2016](https://shiny.rstudio.com/)). By typing 

```r
TSIS.app(data.size.max = 100) 
```
in R console after loading TSIS package, where “data.size.max” is the maximum allowance size for loading input data. The default is 100MB. The App is opened in the default web browser. Users can upload input datasets, set parameters for switch analysis, visualize and save the results easily. The [TSIS](https://github.com/wyguo/TSIS) App includes three tab panels (see <a href="# whole1">Figure 2(A)</a>).

## Tab panel 1: Manual
The first tab panel includes this user manual.

## Tab panel 2: Isoform switch analysis
There are four sections in this panel (see <a href="#fig2">Figure 2</a>).), namely Input data files, Parameter settings, Density/Frequency of switch and output metrics table of isoform switch.


<h2 id="fig2"> </h2>

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_002.png)

**Figure 2:** Second tab panel in TSIS Shiny App. (A) is the three tab panels of the app; (B) is the data input interface; (C) is the interface for TSIS parameter setting; (D) provides the density/frequency plots of isoform switch time and (E) shows the output of TSIS analysis.

### Input data files
Three *.csv format input files can be provided for [TSIS](https://github.com/wyguo/TSIS) analysis. 

- Time-series isoform expression data with first row indicating the replicate labels and second row indicating the time-points. The remaining lines are isoform names in the first column followed by the expression values (see <a href="#fig3">Figure 3(A)</a>).
- Gene and isoform mapping table with gene names in first column  and transcript isoform names in the second column (see <a href="#fig3">Figure 3(B)</a>).
- Optional. A list of isoform names of interest. Users can output subsets of results by limiting the output to a list of isoforms of interest, for example, protein coding transcripts (see <a href="#fig3">Figure 3(C)</a>).


<h2 id="fig3"> </h2>

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_003.png)

**Figure 3:** The format of input csv files for (A) transcript isoform expression, (B) two column table of gene-isoform mapping and (C) A list of isoform names of interest.

<a href="#fig2">Figure 2(B)</a>  and <a href="#fig4">Figure 4(A)</a> shows the data input interface for time-series isoform expression and gene-isoform mapping. By clicking the "Browse…" button, a window is open for data loading (see <a href="#fig4">Figure 4(B)</a>). Users can use the interface shown in <a href="#fig4">Figure 4(C)</a>  to load the names of subset of isoforms. 


<h2 id="fig4"> </h2>

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_004.png)

**Figure 4:** Interface for input information. (A) Input transcript isoform expression and gene-isoform mapping data, (B) is an opened window to select files after clicking “Browser” and (C) is the interface to load isoform names of interest.

### Parameter settings

#### Scoring parameters

The section in <a href="#fig2">Figure 2(C)</a> and <a href="#fig5">Figure 5</a> is used to set the parameters for [TSIS](https://github.com/wyguo/TSIS). The parameters can be set by selecting or typing in corresponding boxes. The details of how to set the parameters are in the text below the “Scoring” button (<a href="#fig5">Figure 5(A)</a>). Scoring data is generated by clicking on the "Scoring" button. Processing tracking bars (<a href="#fig5">Figure 5(B)</a>) are presented at the bottom of the browser when the scoring is in progress.

#### Filtering parameters
<a href="#fig5">Figure 5(C)</a> is the interface for output filtering. Users can set cut-offs, such as for the probability/frequency of switch and sum of average differences, to further refine the switch results. The parameter setting details are in the text under the "Filtering" button.


<h2 id="fig5"> </h2>

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_005.png)

**Figure 5:** TSIS parameter setting section. (A) is the input interface for setting parameters for scoring ; (B) is the process tracking bars and (C) is the interface for setting parameters for filtering.

### Density of switch points
The isoform switches may occur at different time-points in the time-series. To visualize the frequency and density plot of timing of switches, TSIS Shiny App provides the plot interface as shown in <a href="#fig6">Figure 6</a>. Frequency and density bar plots as well as line plots, which correspond to isoform switch time-points after scoring and filtering processes, are presented by clicking the corresponding radio buttons. The plot can be saved in html, pdf or png format.

Note: The plot is made by using [plotly]( https://plot.ly/r/) R package. Users can move the mouse around the plot to show plot values and select part of the plot to zoom in. More actions are available by using the tool bar in the top right corner of the plot.


<h2 id="fig6"> </h2>

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_006.png)

**Figure 6:** Switch time (A) frequency and (B) density plot interface.

### Output results of isoform switch
The output of TSIS analysis can be displayed and exported after scoring or filtering. The columns include the information of isoform names, isoform ratios to genes, the intervals before and after switch, the coordinates of switch points and five measurements to characterize the isoform switch. Table columns can be sorted by clicking the small triangles beside the column names and contents can be searched by typing text in the search box. The explanations for each column are on the top of the table (see <a href="#fig7">Figure 7</a>).


<h2 id="fig7"> </h2>

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_007.png)


**Figure 7:** The output of TSIS.

## Tab panel 3: visualization


<h2 id="fig8"> </h2>

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_008.png)

**Figure 8:** The third tab panel of TSIS Shiny App. (A) is the switch plot section by providing a pair of isoform names. (B) is used to save top n plot into a local folder.

### Switch plots
Any pair of switched transcript isoforms can be visualized by providing their names. Plot type options are error bar plot and ribbon plot (see functions geom_errorbar and geom_smooth in [ggplot2](http://ggplot2.org) package for details) as shown in <a href="#fig8">Figure 8(A)</a> and example plots of G30 in <a href="#fig9">Figure 9</a> and <a href="#fig10">Figure 10</a>. An option is provided to only label the features of switch points with probability/frequency of switch>cut-off in the time frame of interest. The plots can be saved in html ([plotly](https://plot.ly/) format plot), png or pdf format.

### Switch plots in batch
Transcript isoform switch profiles can be plotted in batch by selecting top n (ranking with Score 1 probability/frequency of switch) pairs of isoforms into png or pdf format plots (see <a href="#fig8">Figure 8(B)</a>).

# Running TSIS using command lines -- step by step analysis
In addition to the Shiny App, users can use scripts to do TSIS analysis in R console. The following examples show a step-by-step tutorial of the analysis. Please refer to the function details using help function, e.g. help(iso.switch) or ?iso.switch.

## Loading datasets


```r
##load the data
library(TSIS)
data.exp<-TSIS.data$data.exp
mapping<-TSIS.data$mapping
dim(data.exp);dim(mapping)
```

## Scoring
**Example 1: search intersection points with mean expression**


```r
##Scores
scores.mean2int<-iso.switch(data.exp=data.exp,mapping =mapping,
                     times=rep(1:26,each=9),rank=F,
                     min.t.points =2,min.difference=1,spline =F,
                     spline.df = 9,verbose = F)
```



**Example 2: search intersection points with spline method**

```r
##Scores, set spline=T and define spline degree of freedom to spline.df=9
scores.spline2int<-iso.switch(data.exp=data.exp,mapping =mapping,
                     times=rep(1:26,each=9),rank=F,
                     min.t.points =2,min.difference=1,spline =T,
                     spline.df = 10,verbose = F)
```

## Filtering
**Example 1: general filtering with cut-offs**


```r
##intersection from mean expression
scores.mean2int.filtered<-score.filter(
  scores = scores.mean2int,prob.cutoff = 0.5,diff.cutoff = 1,
  t.points.cutoff = 2,pval.cutoff = 0.001, cor.cutoff = 0,
  data.exp = NULL,mapping = NULL,sub.isoform.list = NULL,
  sub.isoform = F,max.ratio = F,x.value.limit = c(1,26) 
)
scores.mean2int.filtered[1:5,]
```


```r
##intersection from spline method
scores.spline2int.filtered<-score.filter(
  scores = scores.spline2int,prob.cutoff = 0.5,
  diff.cutoff = 1,t.points.cutoff = 2,pval.cutoff = 0.01,
  cor.cutoff = 0.5,data.exp = NULL,mapping = NULL,
  sub.isoform.list = NULL,sub.isoform = F,max.ratio = F,
  x.value.limit = c(9,17) 
)
```



**Example 2: only show subset of results according to an isoform list**


```r
##intersection from mean expression
##input a list of isoform names for investigation.
sub.isoform.list<-TSIS.data$sub.isoforms
sub.isoform.list[1:10]
##assign the isoform name list to sub.isoform.list and set sub.isoform=TRUE
scores.mean2int.filtered.subset<-score.filter(
  scores = scores.mean2int,prob.cutoff = 0.5,diff.cutoff = 1,
  t.points.cutoff = 2,pval.cutoff = 0.01, cor.cutoff = 0.5,
  data.exp = NULL,mapping = NULL,sub.isoform.list = sub.isoform.list,
  sub.isoform = T,max.ratio = F,x.value.limit = c(9,17) 
)
```



**Example 3: only show results of the most abundant transcript within a gene**


```r
scores.mean2int.filtered.maxratio<-score.filter(
  scores = scores.mean2int,prob.cutoff = 0.5,diff.cutoff = 1,
  t.points.cutoff = 2,pval.cutoff = 0.01, cor.cutoff = 0,
  data.exp = data.exp,mapping = mapping,sub.isoform.list = NULL,
  sub.isoform = F,max.ratio = T,x.value.limit = c(9,17) 
)
```

## Make plots

### Switch time points density/frequency


```r
library(gridExtra)
g1<-switch.density(scores.mean2int.filtered$x.value,make.plotly = F,
                   show.line = F,plot.type = 'frequency',
                   title = 'Frequency of switch time' ,time.points = 1:26)
g2<-switch.density(scores.mean2int.filtered$x.value,make.plotly = F,
                   show.line = T,plot.type = 'density',
                   title = 'Density of switch time' ,time.points = 1:26)
gridExtra::grid.arrange(g1,g2,ncol=2)
```

### Error bar plot


<h2 id="fig9"> </h2>

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_009.png)



```r
plotTSIS(data2plot = data.exp,scores = scores.mean2int.filtered,
         iso1 = 'G30.iso2',iso2 = 'G30.iso3',gene.name = NULL,
         y.lab = 'Expression',make.plotly = F,
         times=rep(1:26,each=9),prob.cutoff = 0.5,
         x.lower.boundary = 9,x.upper.boundary = 17,
         show.region = T,show.scores = T,
         line.width =0.5,point.size = 3,
         error.type = 'stderr',show.errorbar = T,errorbar.size = 0.5,
         errorbar.width = 0.2,spline = F,spline.df = NULL,ribbon.plot = F)
```


<h2 id="fig9"> </h2>

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_010.png)


### Ribbon plot


```r
plotTSIS(data2plot = data.exp,scores = scores.mean2int.filtered,
         iso1 = 'G30.iso2',iso2 = 'G30.iso3',gene.name = NULL,
         y.lab = 'Expression',make.plotly = F,
         times=rep(1:26,each=9),prob.cutoff = 0.5,
         x.lower.boundary = 9,x.upper.boundary = 17,
         show.region = T,show.scores = T,
         line.width =0.5,point.size = 3,
         error.type = 'stderr',show.errorbar = T,errorbar.size = 0.5,
         errorbar.width = 0.2,spline = F,spline.df = NULL,ribbon.plot = T)
```


<h2 id="fig10"> </h2>

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/figures_011.png)


# References
Chang, W., et al. 2016. shiny: Web Application Framework for R. https://CRAN.R-project.org/package=shiny

Hastie, T.J. and Tibshirani, R.J. Generalized additive models. Chapter 7 of Statistical Models in S eds. Wadsworth & Brooks/Cole 1992.

Sebestyen, E., Zawisza, M. and Eyras, E. Detection of recurrent alternative splicing switches in tumor samples reveals novel signatures of cancer. Nucleic Acids Res 2015;43(3):1345-1356.





