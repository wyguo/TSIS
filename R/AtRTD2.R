#' Example datasets of AtRTD2.
#'
#' The data sources include subset of TPM time-series transcript expression data with 26 time points, 3 biological replicates and 3
#' technical replicates (234 sample points in total) from the AtRTD2 paper (Zhang, et al., 2016), a gene-isoform mapping table of 300
#' genes (first column) and 766 transcripts (second column) and a subset of 100 isoforms.
#'
#' @docType data
#'
#' @usage data(AtRTD2)
#'
#' @format An object of data list
#'
#' @keywords datasets
#'
#'
#' @references
#' Zhang, R., et al. AtRTD2: A Reference Transcript Dataset for accurate quantification of alternative splicing and expression changes
#' in Arabidopsis thaliana RNA-seq data. bioRxiv 2016.
#'
#' @examples
#' names(AtRTD2)
#' ##expression data
#' AtRTD2$data.exp[1:10,1:5]
#' ##mapping data
#' AtRTD2$mapping[1:10,]
#' ##subset of isoform names
#' AtRTD2$sub.isoforms[1:10]


"AtRTD2"

#' Save AtRTD2 data into csv files
#'
#' @param dir directory to save data. If the directory does not exist, a new folder will be created with the provided name.
#'
#' @return Three csv tables of corresponding data in AtRTD2.
#' @export
#' @examples
#' AtRTD2.example(dir='data')
#'
AtRTD2.example<-function(dir='data'){
  if(!file.exists(dir))
    dir.create(dir)
  sapply(names(AtRTD2),function(x){
    dat<-AtRTD2[[x]]
    if(x=='data.exp'){
      dat<-data.frame(isoforms=rownames(dat),dat,row.names = NULL)
    } else if(x=='sub.isoforms'){
      dat<-data.frame(isoforms=dat,row.names = NULL)
    }
    write.csv(dat,file=paste0(dir,'/',x,'.csv'),row.names = F)

  })
  return('Done!')
}
