#' Isoform switch analysis for time-series data
#'
#' This function is used to search and score isoform switch points in time-series isoform expression data.
#'
#' @details The detailed steps:
#'
#' \bold{Step 1: search for time course intersection points.}
#'
#' The expression for a pair of isoforms \eqn{iso_i} and \eqn{iso_j} may experience a number isoform switch in the whole
#' time duration. Two methods have been included to search for these switch points where the isoforms reverse relative
#' expression profiles.
#' \itemize{
#' \item{\bold{Method 1:}}{ use average expression values across time points. Taking average values of the replicates
#' for time points in the input isoform expression data.}
#' \item{\bold{Method 2:}}{ use nature spline curves to fit the time-series data and find intersection points of the
#' fitted curves for each pair of isoforms. See details in \code{\link{ts.spline}}}
#' }
#'
#' \bold{Step 2: score the isoform switches}
#'
#' We defined 5 parameters to score the quality of isoform switch. The first two are the frequency/probability of switch and the
#' sum of average distance before and after switch, used as Score 1 and Score 2 in iso-kTSP \url{https://bitbucket.org/regulatorygenomicsupf/iso-ktsp}
#' method for two condition comparisons (Sebestyen, et al., 2015).
#' To investigate the switches of two isoforms \eqn{iso_i} and \eqn{iso_j} in two conditions \eqn{c_1} and \eqn{c_2}, \bold{Score 1} is defined as
#' \deqn{S_1(iso_i,iso_j|c_1,c_2)=|p(iso_1>iso2|c_1)+p(iso_1<iso_2|c_2)-1|}
#' where \eqn{p(iso_1>iso2|c_1)} and \eqn{p(iso_1<iso_2|c_2)} are the frequencies/probabilities that the samples of one isoform
#' is greater or less than the other in corresponding conditions. \bold{Score 2} is defined as
#' \deqn{S_2(iso_i,iso_j|c_1,c_2)=|mean.dist(iso_i,iso_2|c_1)|+|mean.dist(ios_1,iso_2|c_2)|}
#' where \eqn{mean.dist(iso_i,iso_2|c_1)} and \eqn{mean.dist(ios_1,iso_2|c_2)} are the mean distances of samples in conditions \eqn{c_1} and \eqn{c_2}, respectively.
#'
#' However, the time-series for a pair of isoforms may undergo a number of switches in the time duration.
#' The time duration is divided into intervals with the intersection points determined in Step 1.
#' To extend the iso-kTSP to TSIS, the samples in each pair of consecutive intervals before and after
#' switch are assimilated as samples in two conditions to implement the calculation of Score 1 and Score 2.
#'
#' The time-series isoform switches are more complex than the comparisons over two conditions. In addition to
#' Score 1 and Score 2 for each switch point, we defined other 3 parameters as metrics of switch qualities.
#'
#' \itemize{
#' \item{\bold{Score 3:}}{ p-value of paired t-test for the two isoform sample differences within each interval.}
#' \item{\bold{Score 4:}}{ Time point number within each interval.}
#' \item{\bold{Score 5:}}{ Pearson correlation of two isoforms.}
#' }
#' Note: since each switch point has a left and right adjoined intervals before and after
#' switch, two p-values and two numbers of time points for the intervals are assigned to each
#' switch point, respectively.
#'
#' @param data.exp isoform expression data frame with row names of isoforms and column names of samples.
#' @param mapping gene-isoform mapping data frame with first column of genes and second column of isoforms.
#' @param t.start,t.end start and end time points of the time-series data. The time step is assumed to be 1.
#' @param nrep number of replicates for each time point.
#' @param rank logical, to use rank of isoform expression for each sample (TRUE) or not (FALSE).
#' @param min.t.points pre-filtering, if the number of time points in all intervals < \code{min.t.points}, skip this pair of isoforms.
#' @param min.distance pre-filtering, if the sample distances in the time courses (mean expression or splined value) for
#' intersection search all < \code{min.distance}, skip this pair of isoforms.
#' @param spline logical, to use spline method (TRUE) or mean expression (FALSE).
#' @param spline.df the degree of freedom used in spline method. See \code{\link{ns}} in \code{\link{splines}} for details.
#' @param verbose logical, to track the progressing of runing (TRUE) or not (FALSE).
#'
#' @references
#' 1.	Sebestyen E, Zawisza M, Eyras E: Detection of recurrent alternative splicing switches in tumor samples reveals novel signatures of cancer.
#' Nucleic Acids Res 2015, 43(3):1345-1356.
#'
#' @return  a data frame of scores. The column names:
#' \itemize{
#' \item{iso1,iso2: }{the isoforms.}
#' \item{iso1.mean.ratio, iso2.mean.ratio: }{the ratio of isoforms to their genes.}
#' \item{left.interval, left.interval: }{the left (before switch) and right (after switch) intervals of a switch point.}
#' \item{x.value, y.value: }{x and y coordinates of switch points.}
#' \item{left.prob, right.prob: }{the frequencies/probabilities that the samples of a isoform is greater or less than the other in left and right intervals, respectively. }
#' \item{left.dist, right.dist: }{the average sample distances in left and right intervals, respectively. }
#' \item{left.pval, right.pval: }{p-values of paired t-test for samples in left and right intervals, respectively. }
#' \item{left.t.points, right.t.points: }{number of time points in left and right intervals, respectively. }
#' \item{prob: }{Score1.}
#' \item{dist: }{Score2.}
#' \item{cor: }{Pearson correlation of two isoforms.}
#'
#' }
#' @export
#'

iso.switch<-function(data.exp,mapping,t.start=1,t.end=26,nrep=9,rank=F,
                     min.t.points=2,min.distance=1,spline=F,spline.df=NULL,verbose=T){


  start.time <- Sys.time()
  ##Step 1: Checking data information
  ##genes more than 2 transcripts
  genes0<-mapping[,1]
  isoforms0<-mapping[,2]

  ##expression ratio aboundance
  data2intersect.ratio<-apply(data.exp,1,mean)
  data2intersect.ratio<-rowratio(x = data2intersect.ratio,group = genes0)

  #choose genes with more than 2 transcripts
  idx<-table(genes0)
  genes<-names(idx)[idx>1]
  isoforms<-isoforms0[which(genes0 %in% genes)]
  data.exp<-data.exp[isoforms,]
  if(rank)
    data2switch<-apply(data.exp,2,rank) else data2switch=data.exp




  message('Input genes: ', length(unique(genes0)))
  message('Genes with more than 2 isoforms: ', length(genes))
  message('Average isoforms per gene for switch analysis: ', round(length(isoforms)/length(genes),3))



  ##step 2: Pickign isoform-pairs with intersection points...
  message(paste0('Step 1: Search for intersection points with ',if(spline) 'Spline method...' else 'Mean expression..'))
  ##Average values of time points
  ##data for searching isoform swtich points
  if(spline){
    if(is.null(spline.df))
      spline.df<-floor(2*(t.end-t.start)/3)
    data2intersect<-t(apply(data.exp,1,function(x) ts.spline(x,t.start = t.start,t.end = t.end,nrep = nrep,df = spline.df)))
  } else {
    data2intersect<-t(rowmean(t(data.exp),group = rep(t.start:t.end,each=nrep)))
  }


  ##Find intersection points
  iso.intersections<-list()
  if(verbose)
    pb <- txtProgressBar(min = 0, max = length(genes), style = 3,width = 75)
  for(i in 1:length(genes)){
    Sys.sleep(0)
    sub.gene<-genes[i]
    sub.isoform<-as.vector(isoforms[which(grepl(sub.gene,isoforms))])

    onegene<-data2intersect[sub.isoform,]
    comb<-combn(sub.isoform,2)
    # iso.intersections<-list()
    for(j in 1:ncol(comb)){
      iso1<-comb[1,j]
      iso2<-comb[2,j]
      ##
      x1=onegene[iso1,]
      x2=onegene[iso2,]

      if(max(abs(x1-x2))<min.distance)
        next

      ##intersection points
      iso.inter<-ts.intersection(x1=x1,x2=x2)
      ##check if have intersection points
      if(is.null(iso.inter))
        next

      ##check if the switching lasting more than 2 consecutive time points
      check.consecutive<-diff(c(0,iso.inter$cross.points$x.points,26))
      if(max(check.consecutive)<min.t.points)
        next

      iso.inter<-data.frame(iso1=iso1,iso2=iso2,t(iso.inter))
      colnames(iso.inter)<-c('iso1','iso2',paste0('switch',1:(ncol(iso.inter)-2)))
      iso.intersections<-c(iso.intersections,setNames(object = list(iso.inter),nm = paste0(iso1,'_vs_',iso2)))
    }
    if(verbose)
      setTxtProgressBar(pb, i)
  }
  if(verbose)
    close(pb)

  message(' ',paste0(length(iso.intersections),' pairs of isoforms have intersection points.'))





  message('Step 2: Calculate scores for isoform switch')
  message(' Score 1: Switch frequencies/probabilities')
  message(' Score 2: Sum of average sample distances before and after switch.')
  message(' Score 3: P-values of sample differences before and after switch')
  message(' Score 4: Time points in each intervals')
  message(' Score 5: Pearson correlation of isoforms')
  ##step 2: calculate probablities

  iso.names<-names(iso.intersections)
  time.points<-rep(t.start:t.end,each=nrep)
  iso.scores<-data.frame()

  if(verbose)
    pb <- txtProgressBar(min = 0, max = length(iso.names), style = 3,width = 75)
  for(i in 1:length(iso.names)){
    Sys.sleep(0)
    #  i=1
    # i=which(grepl('AT1G01060.2_vs_AT1G01060.3',x = iso.names))
    iso<-unlist(strsplit(iso.names[i],'_vs_'))
    iso1<-iso[1]
    iso2<-iso[2]
    #generate dataset for iso for isoform switch analysis
    n.inters<-data.frame(iso.intersections[[iso.names[i]]][,-c(1:2)])
    colnames(n.inters)<-colnames(iso.intersections[[iso.names[i]]])[-c(1:2)]
    inters=as.numeric(n.inters[1,])

    interval=cut(time.points,unique(c(t.start,as.numeric(inters),t.end)),include.lowest = T,dig.lab = 3)
    # sub.data2swith=data.frame(interval=interval,t(data.exp.rank[c(iso1,iso2),]),iso.idx=c(-1,1)[as.numeric(interval)%%2+1])
    sub.data2swith=data.frame(interval=interval,t(data2switch[c(iso1,iso2),]))
    ##add a column of sign of differences

    diff.sign<-sign(sub.data2swith[,2]-sub.data2swith[,3])
    interval.idx<-as.numeric(interval) %%2
    interval.idx[interval.idx==0]<--1
    diff.sign<-diff.sign*interval.idx
    sub.data2swith$diff.sign<-diff.sign
    x0<-diff.sign[diff.sign!=0][1]
    #Define scores

    s<-by(data = sub.data2swith,INDICES =sub.data2swith[,1],FUN = function(x){
      ##score1:
      prob<-length(x[,4][x[,4]==x0])/length(x[,4])
      prob[is.na(prob)]<-0
      ##score2:
      x.diff<-x[,2]-x[,3]
      dist<-mean(x.diff)
      ##score3:
      pval<-t.test(x[,3],x[,2],paired = T)$p.value
      pval[is.na(pval)]<-1
      ##score4:
      inter.length<-nrow(x)/nrep
      x<-data.frame(prob=prob,dist=dist,pval=pval,inter.length=inter.length)
      x
    },simplify = T)
    s<-t(do.call(rbind,s))

    if(ncol(s)<3)
      lf.idx<-1:2 else lf.idx<-c(1,rep(2:(ncol(s)-1),each=2),ncol(s))

    ##score 1: prob
    score1.2side<-matrix(s[1,lf.idx],ncol=2,byrow = T)
    colnames(score1.2side)<-c('left.prob','right.prob')
    score1<-zoo::rollapply(as.numeric(s[1,]), width = 2, FUN = function(x) abs(sum(x)-1))

    ##score 2: distance
    score2<-abs(diff(s[2,]))
    score2.2side<-matrix(s[2,lf.idx],ncol=2,byrow = T)
    colnames(score2.2side)<-c('left.dist','right.dist')

    ##score 3: p-value
    score3<-matrix(s[3,lf.idx],ncol=2,byrow = T)
    colnames(score3)<-c('left.pval','right.pval')

    ##score 4: interval length
    score4<-matrix(s[4,lf.idx],ncol=2,byrow = T)
    colnames(score4)<-c('left.t.points','right.t.points')


    ###average ratio
    iso.mean.ratio<-data.frame(iso1.mean.ratio=data2intersect.ratio[iso1],iso2.mean.ratio=data2intersect.ratio[iso2])

    ###left and right intervals
    inter.lr<-matrix(unique(interval)[lf.idx],ncol=2,byrow = T)
    colnames(inter.lr)<-c('left.interval','right.invertal')


    score<-data.frame(iso1=iso1,iso2=iso2,iso.mean.ratio,inter.lr,x.value=as.numeric(n.inters[1,]),y.value=as.numeric(n.inters[2,]),
                      score1.2side,score2.2side,score3,score4,
                      prob=score1,dist=score2,cor=cor(as.numeric(data.exp[iso1,]),as.numeric(data.exp[iso2,])),row.names = NULL)

    iso.scores<-rbind(score,iso.scores)
    if(verbose)
      setTxtProgressBar(pb, i)
  }
  if(verbose)
    close(pb)
  ##sort
  iso.scores<-iso.scores[order(iso.scores$prob,decreasing = T),]

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message(paste0('Time for analysis: ',round(time.taken,3),' ',attributes(time.taken)$units))

  message('Done!!! ')
  return(data.frame(iso.scores))
}


#' @export
#'

iso.switch.shiny<-function(data.exp,data2intersect=NULL,mapping,t.start=1,t.end=26,nrep=9,rank=F,
                           min.t.points=2,min.distance=1,spline=F,spline.df=5){


  start.time <- Sys.time()
  ##Step 1: Checking data information
  ##genes more than 2 transcripts
  genes0<-mapping[,1]
  isoforms0<-mapping[,2]

  ##expression ratio aboundance
  data2intersect.ratio<-apply(data.exp,1,mean)
  data2intersect.ratio<-rowratio(x = data2intersect.ratio,group = genes0)

  #pick genes with more than 2 transcrits
  idx<-table(genes0)
  genes<-names(idx)[idx>1]
  isoforms<-isoforms0[which(genes0 %in% genes)]
  data.exp<-data.exp[isoforms,]
  if(rank)
    data2switch<-apply(data.exp,2,rank) else data2switch=data.exp




  message('Input genes: ', length(unique(genes0)))
  message('Genes with more than 2 isoforms: ', length(genes))
  message('Average isoforms per gene for switch analysis: ', round(length(isoforms)/length(genes),3))



  ##step 2: Pickign isoform-pairs with intersection points...
  message(paste0('Step 1: Search for intersection points with ',if(spline) 'Spline method...' else 'Mean expression..'))

  ##data for searching isoform swtich points
  if(spline){
    if(is.null(spline.df))
      spline.df<-floor(2*(t.end-t.start)/3)
    data2intersect<-t(apply(data.exp,1,function(x) ts.spline(x,t.start = t.start,t.end = t.end,nrep = nrep,df = spline.df)))
  } else {
    data2intersect<-t(rowmean(t(data.exp),group = rep(t.start:t.end,each=nrep)))
  }



  ##Find intersection points
  iso.intersections<-list()
  withProgress(message = 'Searching intersections: ',value=0,{
    pb <- txtProgressBar(min = 0, max = length(genes), style = 3,width = 75)
    for(i in 1:length(genes)){
      sub.gene<-genes[i]
      sub.isoform<-as.vector(isoforms[which(grepl(sub.gene,isoforms))])

      onegene<-data2intersect[sub.isoform,]
      comb<-combn(sub.isoform,2)
      # iso.intersections<-list()
      for(j in 1:ncol(comb)){
        iso1<-comb[1,j]
        iso2<-comb[2,j]
        ##
        x1=onegene[iso1,]
        x2=onegene[iso2,]

        if(max(abs(x1-x2))<min.distance)
          next

        ##intersection points
        iso.inter<-ts.intersection(x1=x1,x2=x2)
        ##check if have intersection points
        if(is.null(iso.inter))
          next

        ##check if the switching lasting more than 2 consecutive time points
        check.consecutive<-diff(c(0,iso.inter$cross.points$x.points,26))
        if(max(check.consecutive)<min.t.points)
          next

        iso.inter<-data.frame(iso1=iso1,iso2=iso2,t(iso.inter))
        colnames(iso.inter)<-c('iso1','iso2',paste0('switch',1:(ncol(iso.inter)-2)))
        iso.intersections<-c(iso.intersections,setNames(object = list(iso.inter),nm = paste0(iso1,'_vs_',iso2)))
      }

      incProgress(1/length(genes), detail = paste(i, ' of ', length(genes)))
      Sys.sleep(0)
      setTxtProgressBar(pb, i)
    }
    close(pb)
  })

  message(' ',paste0(length(iso.intersections),' pairs of isoforms have intersection points.'))



  message('Step 2: Calculate scores for isoform switch')
  message(' Score 1: Switch frequencies/probabilities')
  message(' Score 2: Sum of average sample distances before and after switch.')
  message(' Score 3: P-values of sample differences before and after switch')
  message(' Score 4: Time points in each intervals')
  message(' Score 5: Pearson correlation of isoforms')
  ##step 2: calculate probablities

  iso.names<-names(iso.intersections)
  time.points<-rep(t.start:t.end,each=nrep)
  iso.scores<-data.frame()
  withProgress(message = 'Scoring isoforms: ',value=0,{
    pb <- txtProgressBar(min = 0, max = length(iso.names), style = 3,width = 75)
    for(i in 1:length(iso.names)){

      #  i=1
      # i=which(grepl('AT1G01060.2_vs_AT1G01060.3',x = iso.names))
      iso<-unlist(strsplit(iso.names[i],'_vs_'))
      iso1<-iso[1]
      iso2<-iso[2]
      #generate dataset for iso for isoform switch analysis
      n.inters<-data.frame(iso.intersections[[iso.names[i]]][,-c(1:2)])
      colnames(n.inters)<-colnames(iso.intersections[[iso.names[i]]])[-c(1:2)]
      inters=as.numeric(n.inters[1,])

      interval=cut(time.points,unique(c(t.start,as.numeric(inters),t.end)),include.lowest = T,dig.lab = 3)
      # sub.data2swith=data.frame(interval=interval,t(data.exp.rank[c(iso1,iso2),]),iso.idx=c(-1,1)[as.numeric(interval)%%2+1])
      sub.data2swith=data.frame(interval=interval,t(data2switch[c(iso1,iso2),]))
      ##add a column of sign of differences

      diff.sign<-sign(sub.data2swith[,2]-sub.data2swith[,3])
      interval.idx<-as.numeric(interval) %%2
      interval.idx[interval.idx==0]<--1
      diff.sign<-diff.sign*interval.idx
      sub.data2swith$diff.sign<-diff.sign
      x0<-diff.sign[diff.sign!=0][1]
      #Define scores

      s<-by(data = sub.data2swith,INDICES =sub.data2swith[,1],FUN = function(x){
        ##score1:
        prob<-length(x[,4][x[,4]==x0])/length(x[,4])
        prob[is.na(prob)]<-0
        ##score2:
        x.diff<-x[,2]-x[,3]
        dist<-mean(x.diff)
        ##score3:
        pval<-t.test(x[,3],x[,2],paired = T)$p.value
        pval[is.na(pval)]<-1
        ##score4:
        inter.length<-nrow(x)/nrep
        x<-data.frame(prob=prob,dist=dist,pval=pval,inter.length=inter.length)
        x
      },simplify = T)
      s<-t(do.call(rbind,s))

      if(ncol(s)<3)
        lf.idx<-1:2 else lf.idx<-c(1,rep(2:(ncol(s)-1),each=2),ncol(s))

      ##score 1: prob
      score1.2side<-matrix(s[1,lf.idx],ncol=2,byrow = T)
      colnames(score1.2side)<-c('left.prob','right.prob')
      score1<-zoo::rollapply(as.numeric(s[1,]), width = 2, FUN = function(x) abs(sum(x)-1))

      ##score 2: distance
      score2<-abs(diff(s[2,]))
      score2.2side<-matrix(s[2,lf.idx],ncol=2,byrow = T)
      colnames(score2.2side)<-c('left.dist','right.dist')

      ##score 3: p-value
      score3<-matrix(s[3,lf.idx],ncol=2,byrow = T)
      colnames(score3)<-c('left.pval','right.pval')

      ##score 4: interval length
      score4<-matrix(s[4,lf.idx],ncol=2,byrow = T)
      colnames(score4)<-c('left.t.points','right.t.points')


      ###average ratio
      iso.mean.ratio<-data.frame(iso1.mean.ratio=data2intersect.ratio[iso1],iso2.mean.ratio=data2intersect.ratio[iso2])

      ###left and right intervals
      inter.lr<-matrix(unique(interval)[lf.idx],ncol=2,byrow = T)
      colnames(inter.lr)<-c('left.interval','right.invertal')


      score<-data.frame(iso1=iso1,iso2=iso2,iso.mean.ratio,inter.lr,x.value=as.numeric(n.inters[1,]),y.value=as.numeric(n.inters[2,]),
                        score1.2side,score2.2side,score3,score4,
                        prob=score1,dist=score2,cor=cor(as.numeric(data.exp[iso1,]),as.numeric(data.exp[iso2,])),row.names = NULL)

      iso.scores<-rbind(score,iso.scores)
      incProgress(1/length(iso.names), detail = paste(i, ' of ', length(iso.names)))
      Sys.sleep(0)
      setTxtProgressBar(pb, i)
    }
    close(pb)
  })
  ##sort
  iso.scores<-iso.scores[order(iso.scores$prob,decreasing = T),]

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message(paste0('Time for analysis: ',round(time.taken,3),' ',attributes(time.taken)$units))

  message('Done!!! ')
  return(data.frame(iso.scores))
}
