#' Filtering the scores for isoform switch
#'
#' Filtering the scores output from \code{\link{iso.switch}}.
#'
#' @details
#' A prospective isoform switch should be:
#'
#' \itemize{
#' \item{}{Have high Score 1 of swtich frequency/probability.}
#' \item{}{With proper value of Score 2 the sum of average distances.}
#' \item{}{The samples in the intervals before and after switch are statistically different.}
#' \item{}{The switch event lasting a few time points in both intervals before and after switch, i.e. the intervals should contain a number of time points.}
#' \item{}{For further details, users can investigate the co-expressed isoform pairs with high Pearson correlation. Note: the isoform pairs with high
#' negative correlation may show better switch pattern if look at the time-series plots.}
#' }
#'
#'  Users may need to investigate subset of isoforms for specific purpose. Three options have been build-in the TSIS package.
#'
#' \itemize{
#' \item{}{Users can set the lower and upper boundaries of a region in the time duration to study the switches only within this region.}
#' \item{}{Users can provide a name list of isoforms to only show the results cantain the isoforms in the list.}
#' \item{}{Users can output subset of results with highest ratios (the proportions of isoforms to the genes) isoforms.}
#' }
#'
#' @param scores the scores output from \code{\link{iso.switch}}.
#' @param prob.cutoff,dist.cutoff,t.points.cutoff,pval.cutoff,cor.cutoff the cut-offs corresponding to switch frequencies/probablities,
#' sum of average distances, p-value and time points cut-offs for both intervals before and after switch and Pearson correlation.
#' @param data.exp,mapping the expression and gene-isoform mapping data.
#' @param sub.isoform.list a vector of isoform names to output subset of the corresponding results.
#' @param sub.isoform logical, to output subset of the results(TRUE) or not (FALSE). If TRUE, \code{sub.isoform.list} must be provided.
#' @param max.ratio logical, to show the subset of results with the isoforms of maximum ratios to the genes (TRUE) or not (FALSE).
#' If TRUE, data.exp and
#' mapping data must be provided to calculate the isoform ratios to the genes using \code{\link{rowratio}}.
#' @param x.value.limit the region of x axis (time) for investigation. If there is no intersection point in this region, the isoform
#' pair is filtered.
#'
#'
#' @export
#'
score.filter<-function(scores,prob.cutoff=0.5,dist.cutoff=1,t.points.cutoff=2,pval.cutoff=0.01,cor.cutoff=0.5,
                       data.exp=NULL,mapping=NULL,sub.isoform.list=NULL,sub.isoform=F,max.ratio=F,
                       x.value.limit=c(9,17)){

  if((is.null(data.exp) | is.null(mapping)) & max.ratio){
    msg<-data.frame('Expression data and mapping data are not provided.',row.names = NULL)
    colnames(msg)<-'Warnings:'
    return(msg)
  } else if(is.null(sub.isoform.list) & sub.isoform){
    msg<-data.frame('Subset of isoform names is not provided.',row.names = NULL)
    colnames(msg)<-'Warnings:'
    return(msg)
  } else {
    scores=scores[which(scores$before.t.points>=t.points.cutoff
                        & scores$after.t.points >= t.points.cutoff
                        & scores$prob>prob.cutoff
                        & scores$dist >dist.cutoff
                        & scores$before.pval <pval.cutoff
                        & scores$after.pval<pval.cutoff
                        & abs(scores$cor)>cor.cutoff
                        & scores$x.value >=x.value.limit[1]
                        & scores$x.value <=x.value.limit[2]),]
    scores<-scores[order(scores[,'prob'],decreasing = T),]
    rownames(scores)<-NULL

    ##subsets of the data

    if(sub.isoform & nrow(scores)>0)
      scores<-scores[which(scores$iso1 %in% sub.isoform.list | scores$iso2 %in% sub.isoform.list),]

    if(max.ratio & nrow(scores)>0){
      genes0<-mapping[,1]
      isoforms0<-mapping[,2]

      ##expression ratio aboundance
      data.exp.mean.ratio<-apply(data.exp,1,mean)
      data.exp.mean.ratio<-data.frame(usage=rowratio(x = data.exp.mean.ratio,group = genes0))
      max.ratio.idx<-by(data.exp.mean.ratio,INDICES = genes0,function(x){
        rownames(x)[x==max(x)]
      },simplify = F)
      max.ratio.idx<-do.call(cbind,max.ratio.idx)
      scores<-scores[which(scores$iso1 %in% max.ratio.idx | scores$iso2 %in% max.ratio.idx),]
    }
    return(scores)

  }


}


