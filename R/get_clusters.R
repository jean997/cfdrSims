
#' Count clusters at several levels of z and set threshold
#'@description Count clusters at several levels of z and set threshold
#'@param stats Matrix of stats produced by one of the get_stats functions.
#'@param pos Length p vector of positions
#'@param bw bandwidth to use for smoothing
#'@param nz Number of theresholds to look at
#'@param level Threshold at which to control FDR.
#' @return A list with elements:
#' \describe{
#' \item{\code{clust_num}}{Matrix of cluster counts nz x (1+n.perm)}
#' \item{\code{bw}}{ Bandwidth used for smoothing}
#' \item{\code{clust}}{ Clusters at chosen level of z }
#' \item{\code{z}}{ All levels looked at}
#' \item{\code{z0}}{ Reference level used for merging}
#' \item{\code{fdr}}{ Estimated fdr at each level}
#' \item{\code{zsel}}{ Selected level of z}
#' }
#'@export
get_clusters <- function(stats, pos, bw, nz, level=0.1){
  N <- ncol(stats)
  x <- ksmooth(pos, stats[,1], bandwidth=bw, x.points=pos)$y
  zmin = quantile(abs(x), probs=0.9)
  z <- seq(zmin, max(abs(x)), length.out=nz)
  z0 = 0.3*zmin

  clust_num = apply(stats, MARGIN=2, FUN=function(st){
	  xs =  ksmooth(pos, st, bandwidth=bw, x.points=pos)$y
    unlist(lapply(z, FUN=function(t){
      count_clusters_merged(x=xs, z=t, z0=z0)$nclust
    }))
	})
  lhat <- rowMeans(clust_num[,2:N])
  fdrhat = lhat/(clust_num[,1]+1)
  if(any(fdrhat < level)){
	  myz = min(z[fdrhat < level])
  }else{
	  myz = Inf
  }
  if(myz < max(z)){
    clust <- count_clusters_merged(x=x, z=myz, z0=z0)
    qI <- clust$clust
  }else{
    qI <- NULL
  }
  R <- list("clust_num" = clust_num,  "bw"=bw, "clust"=qI,
            "z"=z, "z0"=z0, "fdr"=fdrhat,"zsel"=myz)
  return(R)
}
