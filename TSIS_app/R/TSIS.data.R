#' Example datasets of TSIS.data.
#'
#' The TSIS package provides the example datasets "TSIS.data" with 30 genes and 140 isoforms, analyzed in 26 time points, each with 9 replicates. 
#' The isoform expression is in TPM (transcript per million) format. Other type of transcript quantifications, such as read counts, 
#' Percentage Splicing Ins (PSIs) can also be used in TSIS.
#'
#'
#' @docType data
#'
#' @usage data(TSIS.data)
#'
#' @format An object of data list
#'
#' @keywords datasets
#'
#'
#'
#' @examples
#' names(TSIS.data)
#' ##expression data
#' TSIS.data$data.exp[1:10,1:5]
#' ##mapping data
#' TSIS.data$mapping[1:10,]
#' ##subset of isoform names
#' TSIS.data$sub.isoforms[1:10]
#'
#'
 "TSIS.data"


