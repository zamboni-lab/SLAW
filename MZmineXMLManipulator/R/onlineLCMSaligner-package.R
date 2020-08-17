#' Online alignment of LC-MS files
#'
#' Online alignment of LC-MS files
#'
#' Online alignment of LC-MS files
#'
#' @docType package
#' @importFrom ggplot2 ggplot aes geom_point geom_histogram xlab ylab xlim ylim theme
#' @importFrom gghighlight gghighlight
#' @import viridis
#' @import rtree
#' @import Rcpp
#' @import xcms
#' @importFrom BiocParallel bplapply bpmapply bpparam
#' @importFrom stringr str_split fixed
#' @importFrom tools md5sum
#' @importFrom Matrix Matrix
#' @importFrom ClusterR GMM
#' @importFrom lpSolve lp.transport
#' @useDynLib onlineLCMSaligner, .registration = TRUE 
#' @name onlineLCMSaligner-package
NULL