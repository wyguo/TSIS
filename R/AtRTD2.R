#' Example datasets of AtRTD2.
#'
#'The TSIS package provides the example datasets "AtRTD2" with 2,666 genes and 6,307 isoforms,
#'analysed in 26 time points, each with 3 biological replicates and 3 technical replicates.
#'The experiments were designed to investigate the Arabidopsis gene expression response to cold.
#'The isoform expression is in TPM (transcript per million) format. For the experiments and data
#'quantification details, please see the AtRTD2 paper (Zhang, et al.,2016). Other type of
#'transcript quantifications, such as read counts, Percentage Splicing Ins (PSIs) can also be used in TSIS.
#'
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
AtRTD2.example<-function(dir='.'){
  # if(!file.exists(dir))
  #   dir.create(dir)

  download.file('https://github.com/wyguo/TSIS/raw/master/data/example_data.zip',
                destfile = 'example data.zip',quiet = T)
  unzip('example data.zip',exdir = '.')
  invisible(file.remove('example data.zip'))

  # sapply(names(AtRTD2),function(x){
  #   dat<-AtRTD2[[x]]
  #   if(x=='data.exp'){
  #     dat<-data.frame(isoforms=rownames(dat),dat,row.names = NULL)
  #   } else if(x=='sub.isoforms'){
  #     dat<-data.frame(isoforms=dat,row.names = NULL)
  #   }
  #   write.csv(dat,file=paste0(dir,'/',x,'.csv'),row.names = F)
  #
  # })
  return('Done!')
}
