
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
#' \item{\code{clust}}{ Clusters at zsel }
#' \item{\code{z}}{ All levels looked at}
#' \item{\code{z0}}{ Reference level used for merging}
#' \item{\code{fdr}}{ Estimated fdr at each level}
#' \item{\code{zsel}}{ Selected level of z}
#' }
#'@export
get_clusters <- function(stats, pos, bw, nz,
                         level=c(0.02, 0.05, 0.1, 0.2),
                         z.min.quantile=0.9, seg.ends=NULL){
  N <- ncol(stats)
  x <- ksmooth(pos, stats[,1], bandwidth=bw, x.points=pos)$y
  zmin = quantile(abs(x), probs=z.min.quantile)
  z <- seq(zmin, max(abs(x)), length.out=nz)
  z0 <- 0.3*zmin

  clust_num <- apply(stats, MARGIN=2, FUN=function(st){
	  xs =  ksmooth(pos, st, bandwidth=bw, x.points=pos)$y
    count_clusters_merged(x=xs, z=z, z0=z0, seg.ends=seg.ends)
	})
  clust <- list()
  if(is.null(seg.ends)) seg.ends <- c(nrow(stats))
  nseg <- length(seg.ends)
  C <- array(dim=c(nz, nseg, N))
  for(i in 1:N) C[, , i] <- matrix(clust_num[, i], nrow=nz, byrow=FALSE)
  R <- C[, , 1]
  if(nseg==1) R <- matrix(R, nrow=nz, ncol=1)
  lhat <- apply(C[, , 2:N, drop=FALSE], MARGIN=c(1, 2), FUN=mean )
  zsel <- matrix(nrow=nrow(stats), ncol=length(level))
  for(j in 1:length(level)){
    idx <- choose_z(lhat, R, level[j])
    strt <- 1
    for(i in 1:nseg){
      if(!is.finite(idx[i, 1])) zsel[strt:seg.ends[i], j] <- Inf
        else zsel[strt:seg.ends[i], j] <- z[idx[i, 1]]
      strt <- seg.ends[i] + 1
    }

    clust[[j]] <- name_clusters_merged(x=x, z=zsel[,j], z0=z0)

  }
  R <- list("clust_num" = clust_num,  "bw"=bw, "clust"=clust, "x"=x,
            "z"=z, "z0"=z0, "lhat"=lhat,"R"=R, "zsel"=zsel)
  return(R)
}

