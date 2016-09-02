
#' Calculate accuracy rates given a set of regions or a list of statsistics and a threshold
#'@description Calculates tpr, number of true positives, number of false positives, and the false discovery proportion.
#'@param signal Intervals object giving regions of true signal
#'@param discoveries Intervals object giving discoveries
#'@param x vector of statistics. Provide with z and z0 in lieu of discoveries
#'@param z vector of thresholds. Should be length of x
#'@param z0 reference level for merging
#' @return A list with elements:
#' \describe{
#' \item{\code{tpr}}{ True positive rate}
#' \item{\code{nfp}}{ Number of false positives}
#' \item{\code{ntp}}{ Number of true positives}
#' \item{\code{fdp}}{ False discover proportion nfp/(nfp + ntp)}
#' }
#'@export
tpr_nfp <- function(signal, discoveries=NULL, x=NULL, z=NULL, z0=NULL){
  if(is.null(discoveries)){
    if(is.null(x) | is.null(z)) stop("Provide either regions or stats and threshold\n")
    if(is.null(z0)) z0 <- z
    discoveries <- name_clusters_merged(x=x, z=z, z0=z0, labs=FALSE)
  }else if(! "intervals" %in% class(discoveries)){
    discoveries <- Intervals(discoveries)
  }
  if(nrow(signal)==0){
    R <- c("tpr"=NA, "nfp"=nrow(discoveries), "ntp"=0, "fdp"=1)
    return(R)
  }
  if(nrow(discoveries)==0){
    R <- c("tpr"=0, "nfp"=0, "ntp"=0, "fdp"=0)
    return(R)
  }
  if(!"intervals" %in% class(signal)){
    stopifnot(dim(signal)[2]==2)
    signal <- Intervals(signal)
  }
  td <- interval_overlap(discoveries, signal)
  #1 if overlapping with signal, 0 otherwise
  td_class <- unlist(lapply(td, FUN=length))
  td <- unlist(td)

  n.t.disc <- sum(td_class > 0)

  n.f.disc <- sum(td_class==0)


  R <- c("tpr"=length(unique(td))/nrow(signal), "nfp"=n.f.disc,
         "ntp"=n.t.disc, "fdp"=n.f.disc/nrow(discoveries))
  return(R)
}
