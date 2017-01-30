---
title: 'TSIS: an R package to infer time-series isoform switch of alternative splicing'
author: "Wenbin Guo"
date: '`r Sys.Date()`'
output:
  html_document:
    code_folding: show
    fig_caption: yes
    highlight: textmate
    theme: cerulean
    toc: yes
    toc_depth: 4
    toc_float:
      collapsed: false
      smooth_scroll: false
subtitle: User manual
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(cache = TRUE,eval=F,echo = T)
```
#Description
[TSIS](https://github.com/wyguo/TSIS) is an R package for detecting transcript isoform switch for time-series data. Transcript isoform switch occurs when a pair of isoforms reverse the order of expression levels as shown in <a href="#fig1">Figure 1</a>. 

<!-- <img src="fig/manual figures_001.png" width="900px" height="200px" /> -->
<!--![](fig/manual figures_001.png)-->
 
 ![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_001.png) 


**Figure 1:** Isoform switch analysis methods. Expression data with 3 replicates for each condition/time point is simulated for isoforms $iso_i$ and $iso_j$. The points in the plots represent the samples and the black lines are the average of samples. (A) is the iso-kTSP algorithm for comparisons of two conditions $c_1$ and $c_2$. The iso-kTSP is extended to time-series isoform switch (TSIS) in figure (B). The time-series with 6 time points is divided into 4 intervals by the intersection points of average expression. Five features for switch evaluation are determined based on the intervals before and after switch, e.g. the before and after intervals adjoined to switch point $P_i$. <h2 id="fig1"> </h2>
TSIS characterizes the transcript switch by 1) defining the isoform switch time points for any pair of transcript isoforms within a gene, 2) describing the switch using 5 different features, 3) filtering the results with user’s specifications and 4) visualizing the results using different plots for the user to examine further details of the switches. All the functions are available in the forms of a graphic interface implemented by [Shiny App](https://shiny.rstudio.com/) (a web application framework for R) ([Chang, et al., 2016](https://shiny.rstudio.com/)), in which users can implement the analysis as easy as mouse click. The tool can also be run just in command lines without graphic interface. This tutorial will cover both in the following sections. 

##Determine switch points
Given that a pair of isoforms $iso_i$ and $iso_j$ may have a number of switches in a time-series, we have offered two approaches to search in TSIS:

- The first approach takes the average values of the replicates for each time points for each isoform expression data and search for the cross points for the average value of two isoforms across the time points (as the example in <a href="#fig1">Figure1 (B)</a>).
- The second approach uses natural spline curves to fit the time-series data for each transcript isoform and find cross points of the fitted curves for each pair of isoforms. Please see the ns() function in [splines](https://stat.ethz.ch/R-manual/R-devel/library/splines/html/ns.html) package for details of spline fitting of time-series.  

##Define switch scoring features
After switch points are determined, we define each switch by 1) the switch point $P_i$ , 2) time points between switch points $P_{i-1}$  and $P_i$ as interval before switch  $P_i$ and 3) time points between switch points $P_i$ and $P_{i+1}$  as interval after the switch $P_i$ (see <a href="#fig1">Figure 1(B)</a>)
We defined 5 features to score each isoform switch. The first two are the probability/frequency of switch and the sum of average distance before and after switch, used as Score 1 and Score 2 in [iso-kTSP](https://bitbucket.org/regulatorygenomicsupf/iso-ktsp) method [(Sebestyen, et al., 2015)]( http://nar.oxfordjournals.org/content/early/2015/01/10/nar.gku1392.full.pdf) (see <a href="#fig1">Figure 1(A)</a>)). 

- For a switch point $P_i$ of two isoforms $iso_i$ and $iso_j$ with before interval $I_1$ and after interval $I_2$ (see example in <a href="#fig1">Figure1 (B)</a>), Feature 1 is defined as
$$F_1 (iso_i,iso_j |I_1,I_2)=|p(iso_i>iso_j |I_1)+p(iso_i<iso_j |I_2)-1|,$$ 
Where  $p(iso_i>iso_j│c_1 )$ and $p(iso_i<iso_j│c_2 )$ are the frequencies/probabilities that the samples of one isoform is greater or less than in the other in corresponding intervals. 
- Feature 2 is defined as
$$F_2 (iso_i,iso_j |I_1,I_2)=|mean.dist(iso_i,iso_j |I_1)|+|mean.dist(ios_i,iso_j |I_2)|,$$
where $mean.dist(iso_i,iso_j│I_1 )$ and $mean.dist(ios_i,iso_j│I_2 )$ are the mean distances of samples in intervals $I_1$ and $I_2$, respectively.
- Feature 3 is the p-value of paired t-test for the two isoform sample differences within each interval. The dependency R function for testing is t.test(), i.e.
$$F_3 (iso_i,iso_j |I_k )=pval⇐t.test(x=iso_i  \text{ samples in } I_k,y=iso_j  \text{ samples in } I_k,paired=TRUE)$$
Where $k=1,2$ represent the indices of the intervals before and after switch.

- Feature 4 indicates the time points numbers in intervals $I_1$ and $I_2$ . For example the intervals before and after switch point $P_i$ in <a href="#fig1">Figure1 (B)</a> both have 2 time points, i.e.
$$F_4 (iso_i,iso_j│I_k )=2,k=1,2$$
- Feature 5 is the Pearson correlation of two isoforms, i.e. 
$$F_5 (iso_i,iso_j )=cor(\text{samples of } iso_1,\text{samples of  } iso_2 )$$

#Installation and loading

## Check before installation 
Due to an issue with [devtools](https://cran.r-project.org/web/packages/devtools/index.html), if  R software is installed in a directory whose name has space character in it, e.g. in "C:\\Program Files", users may get error message "'C:\\Program' is not recognized as an internal or external command". This issue has to be solved by making sure that R is installed in a directory whose name has no space characters. 
Users can check the R installation location by typing
```{r}
R.home()
```

## Install dependency packages
install.packages(c("shiny", "shinythemes","ggplot2","plotly","zoo","gtools","devtools"), dependencies=TRUE)

## Install TSIS package
Install [iso-kTSP](https://bitbucket.org/regulatorygenomicsupf/iso-ktsp) package from Github using [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package.
```{r}
library(devtools)
devtools::install_github("wyguo/TSIS")
```

## Loading
Once installed, TSIS package can be loaded as normal
```{r}
library(TSIS)
```

#Shiny App -- as easy as mouse click
To make the implement more user friendly, TSIS analysis is integrated into a [Shiny App](https://shiny.rstudio.com/) (a web application framework for R) ([Chang, et al., 2016](https://shiny.rstudio.com/)). By typing 
```{r}
TSIS.app() 
```
in R console after loading TSIS package, the App is opened in the default web browser. Users can upload input datasets, set parameters for switch analysis, visualize and save the results as easy as mouse click. The Shiny App includes three tab panels (see <a href=”#fig2”>Figure 2</a>).
<br>
<br>
<!-- <img src="fig/manual figures_002.png" width="600px" height="50px" /> -->
<!--![](fig/manual figures_002.png)-->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_002.png) 

**Figure 2:** The tab panels in Shiny App. <h2 id="fig2"> </h2>

## Tab panel 1: Manual
The first tab panel includes this user manual.

##Tab panel 2: Isoform switch analysis
There are four sections in this panel.

###Input data files
Three types of information are required for [TSIS](https://github.com/wyguo/TSIS) analysis.

- Information 1: Time-series isoform expression data with each row represents an isoform and each column represents a sample (sample replicates followed by one column to another).
- Information 2: Gene and isoform mapping table corresponding to Dataset 1, with first column of gene names and second column of isoform names.
- Information 3: Optional. Names of subset of isoforms. Users can output subset of results by providing a list of isoforms, for example the protein coding transcripts.

<a href="#fig3">Figure 3 (A)</a> shows the data input interface for time-series isoform expression and gene-isoform mapping. By clicking the  “Browse…” button, a window is open for data loading (<a href=”#fig3”>Figure 3(B)</a>). Users can use the interface shown in <a href=“#fig3”>Figure 3(C)</a>  to load the names of subset of isoforms. Please see the following data examples for data format details. 

<br>
<!-- <img src="fig/manual figures_003.png" width="800px" height="400px" /> -->
<!--![](fig/manual figures_003.png)-->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_003.png) 

**Figure 3: ** Interface for input information.  <h2 id="fig3"> </h2>

The [TSIS](https://github.com/wyguo/TSIS) package provides the example datasets "AtRTD2" with 300 genes and 766 isoforms, analysed in 26 time points, each with 3 biological replicates and 3 technical replicates. The experiments were designed to investigate the Arabidopsis gene expression response to cold. The isoform expression is in TPM (transcript per million) format. For the experiments and data quantification details, please see the AtRTD2 paper  [(Zhang, et al.,2016)](http://biorxiv.org/content/early/2016/05/06/051938). Typing the following command to see data information.
```{r}
##26 time points, 3 biological replicates and 3 technical replicates, in total 234 sample points. 
library(TSIS)
AtRTD2$data.exp[1:10,1:3]
AtRTD2$mapping[1:10,]
AtRTD2$sub.isoforms[1:10]
```
Note: The data loaded into the Shiny App must be in *.csv format for loading convenience. Users can download the [example datasets]( https://github.com/wyguo/TSIS/blob/master/data/example_data.zip) from https://github.com/wyguo/TSIS/tree/master/data or by typing the following codes:
```{r}
AtRTD2.example()
```
The data will be saved in a folder "example data" in the working directory. <a href=“#fig4">Figure 4</a> shows the examples of input data in csv format. 

<!-- <img src="fig/manual figures_004.png" width="600px" height="300px" /> -->
<!--![](fig/manual figures_004.png)-->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_004.png) 

**Figure 4:** the format of input csv files for (A) gene expression, (B) gene-isoform mapping and (C) a subset of isoform names. <h2 id=”fig4”></h2>

###Parameter settings

#### Scoring parameters

This section is used to set the parameters for [TSIS](https://github.com/wyguo/TSIS) features. The parameters can be select or typing in corresponding boxes. Scoring process is starting by clicking the “Scoring” button. The parameter setting details are in the text followed the scoring button. Processing tacking bars for time-series intersection points searching (<a href=“#fig3”>Figure 3(B)</a>) and switch scoring (<a href=“#fig3”>Figure 3(C)</a>) for the isoform pairs will present in the bottom of the browser. 

<br>
<!-- <img src="fig/manual figures_005.png" width="800px" height="400px" /> -->
<!--![](fig/manual figures_005.png)-->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_005.png) 

**Figure 5:** Scoring parameter input interface and the processing tracking bars. <h2 id=”fig5”></h2>

#### Filtering parameters
<a href=“#fig6”>Figure 6</a> is the interface for scoring feature filtering. Users can set cut-offs, such as for the probability/frequency of switch and sum of average distances, to further refine the switch results. The parameter setting details are in the text under the “Filtering” button. 

<br>
<!-- <img src="fig/manual figures_006.png" width="800px" height="400px" /> -->
<!--![](fig/manual figures_006.png)-->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_006.png) 

**Figure 6:** the interface for isoform switch. <h2 id=”fig6”></h2>

## Density of switch points
The isoform switches occur at different time points in the time-series. To visualize the frequency and density plot of switch time, TSIS Shiny App provides the plot interface as shown in  <a href=“#fig7”>Figure 7</a>.  Frequency and density bar plots and line plots, which correspond to the “x.value” of switch time column in the following output table, will present by clicking the corresponding radio buttons.  The plot can be saved in html, pdf and png format.

Note: The plot is made by using [plotly]( https://plot.ly/r/) R package. Users can move the mouse around the plot to show plot values and select part of the plot to zoom out. More actions are available by using the tool bar in the top right corner of the plot.

<!-- <img src="fig/manual figures_007.png" width="800px" height="400px" /> -->
<!--![](fig/manual figures_007.png)-->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_007.png) 

**Figure 7:** Switch time density and frequency plot interface. 

## Output scores of isoform switch
The output table of switch features after scoring or filtering. The columns include the information of isoform names, isoform ratios to genes, the intervals before and after switch, the coordinates of switch points and five features of switch quality. Table columns can be sorted by clicking the small triangles beside the column names and contexts can be searched by typing text in the search box. The explanations for each column are on the top of the table (see <a href=“#fig8”>Figure 8</a>). 

<!-- <img src="fig/manual figures_008.png" width="800px" height="400px" /> -->
<!-- ![](fig/manual figures_008.png) -->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_008.png) 

**Figure 8:** The output score feature table. 

#Tab panel 3: Switch visualization
##Switch plots
This part is used to make a time-series plot of a pair of isoforms by providing their names. Plot type options are error bar plot and ribbon plot as shown in <a href=”#fig9”>Figure 9</a> and example plots of <a href=”# AT5G60930_P2_1”> AT5G60930</a> (see functions geom_errorbar and geom_smooth in [ggplot2](http://ggplot2.org) package for details) package for details). An option is provided to only label the features of switch points with probability/frequency of switch>cut-off in the time region for investigation. The plots can be saved in html ([plotly](https://plot.ly/) format plot), png or pdf format (see examples in <a href=“#fig9”>Figure 9</a>, Error bar plot and Ribbon plot).

<!-- <img src="fig/manual figures_009.png" width="800px" height="400px" /> -->
<!-- ![](fig/manual figures_009.png)-->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_009.png) 

**Figure 9**: Isoform switch plot interface. <h2 id=”fig9”><h2>
##Multiple plots for switch
This section is used to save top n (ranking with Feature 1 probability/frequency of switch) pairs of isoforms into png or pdf format plots (see <a href=”#fig10”>Figure 10</a>).

<!-- <img src="fig/manual figures_010.png" width="800px" height="400px" /> -->
<!-- ![](fig/manual figures_010.png) -->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_010.png)

**Figure 9**: Multiple isoform switch plot interface. <h2 id=”fig9”><h2>

#TSIS scripts -- step by step analysis
In addition to the Shiny App, users can use scripts to do TSIS analysis in R console. The following examples show a step-by-step tutorial of the analysis. Please refer to the function details using help function, e.g. help(iso.switch) or ?iso.switch.

## Loading datasets

```{r,eval=T,echo=T}
##load the data
library(TSIS)
data.exp<-AtRTD2$data.exp
mapping<-AtRTD2$mapping
dim(data.exp);dim(mapping)
```

## Scoring
**Example 1: search intersection points with mean expression**

```{r,echo=T,eval=T}
##Scores
scores.mean2int<-iso.switch(data.exp=data.exp,mapping =mapping,
                     t.start=1,t.end=26,nrep=9,rank=F,
                     min.t.points =2,min.distance=1,spline =F,spline.df = 9,verbose = F)
```
<br>
**Example 2: search intersection points with spline method**
```{r}
##Scores, set spline=T and define spline degree of freedom to spline.df=9scores.spline2int<-iso.switch(data.exp=data.exp,mapping =mapping,
                     t.start=1,t.end=26,nrep=9,rank=F,
                     min.t.points =2,min.distance=1,spline =T,spline.df = 9,verbose = F)

```

##Filtering
**Example 1, general filtering with cut-offs**

```{r,eval=T}
##intersection from mean expression
scores.mean2int.filtered<-score.filter(scores = scores.mean2int,prob.cutoff = 0.5,dist.cutoff = 1,
                                       t.points.cutoff = 2,pval.cutoff = 0.01, cor.cutoff = 0.5,
                                       data.exp = NULL,mapping = NULL,sub.isoform.list = NULL,
                                       sub.isoform = F,max.ratio = F,x.value.limit = c(9,17) )

scores.mean2int.filtered[1:5,]
```

```{r}
##intersection from spline method
scores.spline2int.filtered<-score.filter(scores = scores.spline2int,prob.cutoff = 0.5,
                                         dist.cutoff = 1,t.points.cutoff = 2,pval.cutoff = 0.01,
                                         cor.cutoff = 0.5,data.exp = NULL,mapping = NULL,
                                         sub.isoform.list = NULL,sub.isoform = F,max.ratio = F,
                                         x.value.limit = c(9,17) )
```

<br>
 **Example 2, only show subset of results according to an isoform list**

```{r)
##intersection from mean expression
##input a list of isoform names for investigation.
sub.isoform.list<-AtRTD2$sub.isoforms
sub.isoform.list[1:10]
##assign the isoform name list to sub.isoform.list and set sub.isoform=TRUE
scores.mean2int.filtered.subset<-score.filter(scores = scores.mean2int,prob.cutoff = 0.5,dist.cutoff = 1,
                                       t.points.cutoff = 2,pval.cutoff = 0.01, cor.cutoff = 0.5,
                                       data.exp = NULL,mapping = NULL,sub.isoform.list = sub.isoform.list,
                                       sub.isoform = T,max.ratio = F,x.value.limit = c(9,17) )
```

<br>
 **Example 3, only show results of isoforms of maximum ratios to genes**

```{r}
scores.mean2int.filtered.maxratio<-score.filter(scores = scores.mean2int,prob.cutoff = 0.5,dist.cutoff = 1,
                                       t.points.cutoff = 2,pval.cutoff = 0.01, cor.cutoff = 0,
                                       data.exp = data.exp,mapping = mapping,sub.isoform.list = NULL,
                                       sub.isoform = F,max.ratio = T,x.value.limit = c(9,17) )
```

##Make plots

###Error bar plot

```{r}
plotTSIS(data2plot = data.exp,scores = scores.mean2int.filtered,iso1 = 'AT5G60930_P2',
        iso2 = ' AT5G60930_P3',gene.name = NULL,y.lab = 'Expression',make.plotly = F,
        t.start = 1,t.end = 26,nrep = 9,prob.cutoff = 0.5,x.lower.boundary = 9,
        x.upper.boundary = 17,show.region = T,show.scores = T,
        line.width =0.5,point.size = 3,error.type = 'stderr',show.errorbar = T,errorbar.size = 0.5,
        errorbar.width = 0.2,spline = F,spline.df = NULL,ribbon.plot = F )
```


<!-- <img src="fig/AT5G60930_P2_1.png" width="800px" height="400px" /> -->
<!-- ![](fig/AT5G60930_P2_1.png) -->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/AT5G60930_P2_1.png)
<h2 id=”AT5G60930_P2_1”><h2>

###Ribbon plot

```{r}
plotTSIS(data2plot = data.exp,scores = scores.mean2int.filtered,iso1 = 'AT5G60930_P2',
        iso2 = 'AT5G60930_P3',gene.name = NULL,y.lab = 'Expression',make.plotly = F,
        t.start = 1,t.end = 26,nrep = 9,prob.cutoff = 0.5,x.lower.boundary = 9,
        x.upper.boundary = 17,show.region = T,show.scores = T,error.type = 'stderr',
        line.width =0.5,point.size = 3,show.errorbar = T,errorbar.size = 0.5,
        errorbar.width = 0.2,spline = F,spline.df = NULL,ribbon.plot = T )

```

<!-- <img src="fig/AT5G60930_P2_2.png" width="800px" height="400px" /> -->
<!-- ![](fig/AT5G60930_P2_2.png) -->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/AT5G60930_P2_2.png)
<h2 id=”AT5G60930_P2_2”><h2>

#References
Chang, W., et al. 2016. shiny: Web Application Framework for R. https://CRAN.R-project.org/package=shiny

Sebestyen, E., Zawisza, M. and Eyras, E. Detection of recurrent alternative splicing switches in tumor samples reveals novel signatures of cancer. Nucleic Acids Res 2015;43(3):1345-1356.

Zhang, R., et al. AtRTD2: A Reference Transcript Dataset for accurate quantification of alternative splicing and expression changes in Arabidopsis thaliana RNA-seq data. bioRxiv 2016.


=======
---
title: 'TSIS: an R package to infer time-series isoform switch of alternative splicing'
author: "Wenbin Guo"
date: '`r Sys.Date()`'
output:
  html_document:
    code_folding: show
    fig_caption: yes
    highlight: textmate
    theme: cerulean
    toc: yes
    toc_depth: 4
    toc_float:
      collapsed: false
      smooth_scroll: false
subtitle: User manual
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(cache = TRUE,eval=F,echo = T)
```
#Description
[TSIS](https://github.com/wyguo/TSIS) is an R package for detecting transcript isoform switch for time-series data. Transcript isoform switch occurs when a pair of isoforms reverse the order of expression levels as shown in <a href="#fig1">Figure 1</a>. 

<!-- <img src="fig/manual figures_001.png" width="900px" height="200px" /> -->
<!--![](fig/manual figures_001.png)-->
 
 ![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_001.png) 


**Figure 1:** Isoform switch analysis methods. Expression data with 3 replicates for each condition/time point is simulated for isoforms $iso_i$ and $iso_j$. The points in the plots represent the samples and the black lines are the average of samples. (A) is the iso-kTSP algorithm for comparisons of two conditions $c_1$ and $c_2$. The iso-kTSP is extended to time-series isoform switch (TSIS) in figure (B). The time-series with 6 time points is divided into 4 intervals by the intersection points of average expression. Five features for switch evaluation are determined based on the intervals before and after switch, e.g. the before and after intervals adjoined to switch point $P_i$. <h2 id="fig1"> </h2>
TSIS characterizes the transcript switch by 1) defining the isoform switch time points for any pair of transcript isoforms within a gene, 2) describing the switch using 5 different features, 3) filtering the results with user’s specifications and 4) visualizing the results using different plots for the user to examine further details of the switches. All the functions are available in the forms of a graphic interface implemented by [Shiny App](https://shiny.rstudio.com/) (a web application framework for R) ([Chang, et al., 2016](https://shiny.rstudio.com/)), in which users can implement the analysis as easy as mouse click. The tool can also be run just in command lines without graphic interface. This tutorial will cover both in the following sections. 

##Determine switch points
Given that a pair of isoforms $iso_i$ and $iso_j$ may have a number of switches in a time-series, we have offered two approaches to search in TSIS:

- The first approach takes the average values of the replicates for each time points for each isoform expression data and search for the cross points for the average value of two isoforms across the time points (as the example in <a href="#fig1">Figure1 (B)</a>).
- The second approach uses natural spline curves to fit the time-series data for each transcript isoform and find cross points of the fitted curves for each pair of isoforms. Please see the ns() function in [splines](https://stat.ethz.ch/R-manual/R-devel/library/splines/html/ns.html) package for details of spline fitting of time-series.  

##Define switch scoring features
After switch points are determined, we define each switch by 1) the switch point $P_i$ , 2) time points between switch points $P_{i-1}$  and $P_i$ as interval before switch  $P_i$ and 3) time points between switch points $P_i$ and $P_{i+1}$  as interval after the switch $P_i$ (see <a href="#fig1">Figure 1(B)</a>)
We defined 5 features to score each isoform switch. The first two are the probability/frequency of switch and the sum of average distance before and after switch, used as Score 1 and Score 2 in [iso-kTSP](https://bitbucket.org/regulatorygenomicsupf/iso-ktsp) method [(Sebestyen, et al., 2015)]( http://nar.oxfordjournals.org/content/early/2015/01/10/nar.gku1392.full.pdf) (see <a href="#fig1">Figure 1(A)</a>)). 

- For a switch point $P_i$ of two isoforms $iso_i$ and $iso_j$ with before interval $I_1$ and after interval $I_2$ (see example in <a href="#fig1">Figure1 (B)</a>), Feature 1 is defined as
$$F_1 (iso_i,iso_j |I_1,I_2)=|p(iso_i>iso_j |I_1)+p(iso_i<iso_j |I_2)-1|,$$ 
Where  $p(iso_i>iso_j│c_1 )$ and $p(iso_i<iso_j│c_2 )$ are the frequencies/probabilities that the samples of one isoform is greater or less than in the other in corresponding intervals. 
- Feature 2 is defined as
$$F_2 (iso_i,iso_j |I_1,I_2)=|mean.dist(iso_i,iso_j |I_1)|+|mean.dist(ios_i,iso_j |I_2)|,$$
where $mean.dist(iso_i,iso_j│I_1 )$ and $mean.dist(ios_i,iso_j│I_2 )$ are the mean distances of samples in intervals $I_1$ and $I_2$, respectively.
- Feature 3 is the p-value of paired t-test for the two isoform sample differences within each interval. The dependency R function for testing is t.test(), i.e.
$$F_3 (iso_i,iso_j |I_k )=pval⇐t.test(x=iso_i  \text{ samples in } I_k,y=iso_j  \text{ samples in } I_k,paired=TRUE)$$
Where $k=1,2$ represent the indices of the intervals before and after switch.

- Feature 4 indicates the time points numbers in intervals $I_1$ and $I_2$ . For example the intervals before and after switch point $P_i$ in <a href="#fig1">Figure1 (B)</a> both have 2 time points, i.e.
$$F_4 (iso_i,iso_j│I_k )=2,k=1,2$$
- Feature 5 is the Pearson correlation of two isoforms, i.e. 
$$F_5 (iso_i,iso_j )=cor(\text{samples of } iso_1,\text{samples of  } iso_2 )$$

#Installation and loading

## Check before installation 
Due to an issue with [devtools](https://cran.r-project.org/web/packages/devtools/index.html), if  R software is installed in a directory whose name has space character in it, e.g. in "C:\\Program Files", users may get error message "'C:\\Program' is not recognized as an internal or external command". This issue has to be solved by making sure that R is installed in a directory whose name has no space characters. 
Users can check the R installation location by typing
```{r}
R.home()
```

## Install dependency packages
install.packages(c("shiny", "shinythemes","ggplot2","plotly","zoo","gtools","devtools"), dependencies=TRUE)

## Install TSIS package
Install [iso-kTSP](https://bitbucket.org/regulatorygenomicsupf/iso-ktsp) package from Github using [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package.
```{r}
library(devtools)
devtools::install_github("wyguo/TSIS")
```

## Loading
Once installed, TSIS package can be loaded as normal
```{r}
library(TSIS)
```

#Shiny App -- as easy as mouse click
To make the implement more user friendly, TSIS analysis is integrated into a [Shiny App](https://shiny.rstudio.com/) (a web application framework for R) ([Chang, et al., 2016](https://shiny.rstudio.com/)). By typing 
```{r}
TSIS.app() 
```
in R console after loading TSIS package, the App is opened in the default web browser. Users can upload input datasets, set parameters for switch analysis, visualize and save the results as easy as mouse click. The Shiny App includes three tab panels (see <a href=”#fig2”>Figure 2</a>).
<br>
<br>
<!-- <img src="fig/manual figures_002.png" width="600px" height="50px" /> -->
<!--![](fig/manual figures_002.png)-->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_002.png) 

**Figure 2:** The tab panels in Shiny App. <h2 id="fig2"> </h2>

## Tab panel 1: Manual
The first tab panel includes this user manual.

##Tab panel 2: Isoform switch analysis
There are four sections in this panel.

###Input data files
Three types of information are required for [TSIS](https://github.com/wyguo/TSIS) analysis.

- Information 1: Time-series isoform expression data with each row represents an isoform and each column represents a sample (sample replicates followed by one column to another).
- Information 2: Gene and isoform mapping table corresponding to Dataset 1, with first column of gene names and second column of isoform names.
- Information 3: Optional. Names of subset of isoforms. Users can output subset of results by providing a list of isoforms, for example the protein coding transcripts.

<a href="#fig3">Figure 3 (A)</a> shows the data input interface for time-series isoform expression and gene-isoform mapping. By clicking the  “Browse…” button, a window is open for data loading (<a href=”#fig3”>Figure 3(B)</a>). Users can use the interface shown in <a href=“#fig3”>Figure 3(C)</a>  to load the names of subset of isoforms. Please see the following data examples for data format details. 

<br>
<!-- <img src="fig/manual figures_003.png" width="800px" height="400px" /> -->
<!--![](fig/manual figures_003.png)-->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_003.png) 

**Figure 3: ** Interface for input information.  <h2 id="fig3"> </h2>

The [TSIS](https://github.com/wyguo/TSIS) package provides the example datasets "AtRTD2" with 300 genes and 766 isoforms, analysed in 26 time points, each with 3 biological replicates and 3 technical replicates. The experiments were designed to investigate the Arabidopsis gene expression response to cold. The isoform expression is in TPM (transcript per million) format. For the experiments and data quantification details, please see the AtRTD2 paper  [(Zhang, et al.,2016)](http://biorxiv.org/content/early/2016/05/06/051938). Typing the following command to see data information.
```{r}
##26 time points, 3 biological replicates and 3 technical replicates, in total 234 sample points. 
library(TSIS)
AtRTD2$data.exp[1:10,1:3]
AtRTD2$mapping[1:10,]
AtRTD2$sub.isoforms[1:10]
```
Note: The data loaded into the Shiny App must be in *.csv format for loading convenience. Users can download the [example datasets]( https://github.com/wyguo/TSIS/blob/master/data/example_data.zip) from https://github.com/wyguo/TSIS/tree/master/data or by typing the following codes:
```{r}
AtRTD2.example()
```
The data will be saved in a folder "example data" in the working directory. <a href=“#fig4">Figure 4</a> shows the examples of input data in csv format. 

<!-- <img src="fig/manual figures_004.png" width="600px" height="300px" /> -->
<!--![](fig/manual figures_004.png)-->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_004.png) 

**Figure 4:** the format of input csv files for (A) gene expression, (B) gene-isoform mapping and (C) a subset of isoform names. <h2 id=”fig4”></h2>

###Parameter settings

#### Scoring parameters

This section is used to set the parameters for [TSIS](https://github.com/wyguo/TSIS) features. The parameters can be select or typing in corresponding boxes. Scoring process is starting by clicking the “Scoring” button. The parameter setting details are in the text followed the scoring button. Processing tacking bars for time-series intersection points searching (<a href=“#fig3”>Figure 3(B)</a>) and switch scoring (<a href=“#fig3”>Figure 3(C)</a>) for the isoform pairs will present in the bottom of the browser. 

<br>
<!-- <img src="fig/manual figures_005.png" width="800px" height="400px" /> -->
<!--![](fig/manual figures_005.png)-->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_005.png) 

**Figure 5:** Scoring parameter input interface and the processing tracking bars. <h2 id=”fig5”></h2>

#### Filtering parameters
<a href=“#fig6”>Figure 6</a> is the interface for scoring feature filtering. Users can set cut-offs, such as for the probability/frequency of switch and sum of average distances, to further refine the switch results. The parameter setting details are in the text under the “Filtering” button. 

<br>
<!-- <img src="fig/manual figures_006.png" width="800px" height="400px" /> -->
<!--![](fig/manual figures_006.png)-->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_006.png) 

**Figure 6:** the interface for isoform switch. <h2 id=”fig6”></h2>

## Density of switch points
The isoform switches occur at different time points in the time-series. To visualize the frequency and density plot of switch time, TSIS Shiny App provides the plot interface as shown in  <a href=“#fig7”>Figure 7</a>.  Frequency and density bar plots and line plots, which correspond to the “x.value” of switch time column in the following output table, will present by clicking the corresponding radio buttons.  The plot can be saved in html, pdf and png format.

Note: The plot is made by using [plotly]( https://plot.ly/r/) R package. Users can move the mouse around the plot to show plot values and select part of the plot to zoom out. More actions are available by using the tool bar in the top right corner of the plot.

<!-- <img src="fig/manual figures_007.png" width="800px" height="400px" /> -->
<!--![](fig/manual figures_007.png)-->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_007.png) 

**Figure 7:** Switch time density and frequency plot interface. 

## Output scores of isoform switch
The output table of switch features after scoring or filtering. The columns include the information of isoform names, isoform ratios to genes, the intervals before and after switch, the coordinates of switch points and five features of switch quality. Table columns can be sorted by clicking the small triangles beside the column names and contexts can be searched by typing text in the search box. The explanations for each column are on the top of the table (see <a href=“#fig8”>Figure 8</a>). 

<!-- <img src="fig/manual figures_008.png" width="800px" height="400px" /> -->
<!-- ![](fig/manual figures_008.png) -->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_008.png) 

**Figure 8:** The output score feature table. 

#Tab panel 3: Switch visualization
##Switch plots
This part is used to make a time-series plot of a pair of isoforms by providing their names. Plot type options are error bar plot and ribbon plot as shown in <a href=”#fig9”>Figure 9</a> and example plots of <a href=”# AT5G60930_P2_1”> AT5G60930</a> (see functions geom_errorbar and geom_smooth in [ggplot2](http://ggplot2.org) package for details) package for details). An option is provided to only label the features of switch points with probability/frequency of switch>cut-off in the time region for investigation. The plots can be saved in html ([plotly](https://plot.ly/) format plot), png or pdf format (see examples in <a href=“#fig9”>Figure 9</a>, Error bar plot and Ribbon plot).

<!-- <img src="fig/manual figures_009.png" width="800px" height="400px" /> -->
<!-- ![](fig/manual figures_009.png)-->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_009.png) 

**Figure 9**: Isoform switch plot interface. <h2 id=”fig9”><h2>
##Multiple plots for switch
This section is used to save top n (ranking with Feature 1 probability/frequency of switch) pairs of isoforms into png or pdf format plots (see <a href=”#fig10”>Figure 10</a>).

<!-- <img src="fig/manual figures_010.png" width="800px" height="400px" /> -->
<!-- ![](fig/manual figures_010.png) -->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/manual figures_010.png)

**Figure 9**: Multiple isoform switch plot interface. <h2 id=”fig9”><h2>

#TSIS scripts -- step by step analysis
In addition to the Shiny App, users can use scripts to do TSIS analysis in R console. The following examples show a step-by-step tutorial of the analysis. Please refer to the function details using help function, e.g. help(iso.switch) or ?iso.switch.

## Loading datasets

```{r,eval=T,echo=T}
##load the data
library(TSIS)
data.exp<-AtRTD2$data.exp
mapping<-AtRTD2$mapping
dim(data.exp);dim(mapping)
```

## Scoring
**Example 1: search intersection points with mean expression**

```{r,echo=T,eval=T}
##Scores
scores.mean2int<-iso.switch(data.exp=data.exp,mapping =mapping,
                     t.start=1,t.end=26,nrep=9,rank=F,
                     min.t.points =2,min.distance=1,spline =F,spline.df = 9,verbose = F)
```
<br>
**Example 2: search intersection points with spline method**
```{r}
##Scores, set spline=T and define spline degree of freedom to spline.df=9scores.spline2int<-iso.switch(data.exp=data.exp,mapping =mapping,
                     t.start=1,t.end=26,nrep=9,rank=F,
                     min.t.points =2,min.distance=1,spline =T,spline.df = 9,verbose = F)

```

##Filtering
**Example 1, general filtering with cut-offs**

```{r,eval=T}
##intersection from mean expression
scores.mean2int.filtered<-score.filter(scores = scores.mean2int,prob.cutoff = 0.5,dist.cutoff = 1,
                                       t.points.cutoff = 2,pval.cutoff = 0.01, cor.cutoff = 0.5,
                                       data.exp = NULL,mapping = NULL,sub.isoform.list = NULL,
                                       sub.isoform = F,max.ratio = F,x.value.limit = c(9,17) )

scores.mean2int.filtered[1:5,]
```

```{r}
##intersection from spline method
scores.spline2int.filtered<-score.filter(scores = scores.spline2int,prob.cutoff = 0.5,
                                         dist.cutoff = 1,t.points.cutoff = 2,pval.cutoff = 0.01,
                                         cor.cutoff = 0.5,data.exp = NULL,mapping = NULL,
                                         sub.isoform.list = NULL,sub.isoform = F,max.ratio = F,
                                         x.value.limit = c(9,17) )
```

<br>
 **Example 2, only show subset of results according to an isoform list**

```{r)
##intersection from mean expression
##input a list of isoform names for investigation.
sub.isoform.list<-AtRTD2$sub.isoforms
sub.isoform.list[1:10]
##assign the isoform name list to sub.isoform.list and set sub.isoform=TRUE
scores.mean2int.filtered.subset<-score.filter(scores = scores.mean2int,prob.cutoff = 0.5,dist.cutoff = 1,
                                       t.points.cutoff = 2,pval.cutoff = 0.01, cor.cutoff = 0.5,
                                       data.exp = NULL,mapping = NULL,sub.isoform.list = sub.isoform.list,
                                       sub.isoform = T,max.ratio = F,x.value.limit = c(9,17) )
```

<br>
 **Example 3, only show results of isoforms of maximum ratios to genes**

```{r}
scores.mean2int.filtered.maxratio<-score.filter(scores = scores.mean2int,prob.cutoff = 0.5,dist.cutoff = 1,
                                       t.points.cutoff = 2,pval.cutoff = 0.01, cor.cutoff = 0,
                                       data.exp = data.exp,mapping = mapping,sub.isoform.list = NULL,
                                       sub.isoform = F,max.ratio = T,x.value.limit = c(9,17) )
```

##Make plots

###Error bar plot

```{r}
plotTSIS(data2plot = data.exp,scores = scores.mean2int.filtered,iso1 = 'AT5G60930_P2',
        iso2 = ' AT5G60930_P3',gene.name = NULL,y.lab = 'Expression',make.plotly = F,
        t.start = 1,t.end = 26,nrep = 9,prob.cutoff = 0.5,x.lower.boundary = 9,
        x.upper.boundary = 17,show.region = T,show.scores = T,
        line.width =0.5,point.size = 3,error.type = 'stderr',show.errorbar = T,errorbar.size = 0.5,
        errorbar.width = 0.2,spline = F,spline.df = NULL,ribbon.plot = F )
```


<!-- <img src="fig/AT5G60930_P2_1.png" width="800px" height="400px" /> -->
<!-- ![](fig/AT5G60930_P2_1.png) -->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/AT5G60930_P2_1.png)
<h2 id=”AT5G60930_P2_1”><h2>
###Ribbon plot

```{r}
plotTSIS(data2plot = data.exp,scores = scores.mean2int.filtered,iso1 = 'AT5G60930_P2',
        iso2 = 'AT5G60930_P3',gene.name = NULL,y.lab = 'Expression',make.plotly = F,
        t.start = 1,t.end = 26,nrep = 9,prob.cutoff = 0.5,x.lower.boundary = 9,
        x.upper.boundary = 17,show.region = T,show.scores = T,error.type = 'stderr',
        line.width =0.5,point.size = 3,show.errorbar = T,errorbar.size = 0.5,
        errorbar.width = 0.2,spline = F,spline.df = NULL,ribbon.plot = T )

```

<!-- <img src="fig/AT5G60930_P2_2.png" width="800px" height="400px" /> -->
<!-- ![](fig/AT5G60930_P2_2.png) -->

![](https://github.com/wyguo/TSIS/blob/master/vignettes/fig/AT5G60930_P2_2.png)
<h2 id=”AT5G60930_P2_2”><h2>

#References
Chang, W., et al. 2016. shiny: Web Application Framework for R. https://CRAN.R-project.org/package=shiny

Sebestyen, E., Zawisza, M. and Eyras, E. Detection of recurrent alternative splicing switches in tumor samples reveals novel signatures of cancer. Nucleic Acids Res 2015;43(3):1345-1356.

Zhang, R., et al. AtRTD2: A Reference Transcript Dataset for accurate quantification of alternative splicing and expression changes in Arabidopsis thaliana RNA-seq data. bioRxiv 2016.


>>>>>>> origin/master
