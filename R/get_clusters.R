
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
get_clusters <- function(stats, pos, bw, nlam, lambda.max=NULL,
                         signal=NULL, level=c(0.02, 0.05, 0.1, 0.2),
                         z.min.quantile=0.9, seg.ends=NULL){
  N <- ncol(stats)
  x <- ksmooth(pos, stats[,1], bandwidth=bw, x.points=pos)$y
  zmin = quantile(abs(x), probs=z.min.quantile)
  z0 <- 0.3*zmin
  p <- dim(stats)[1]
  if(is.null(seg.ends)) seg.ends <- c(p)
  nseg <- length(seg.ends)
  d <- c(seg.ends[1], diff(seg.ends))

  z <- choose_z_even(stats[, 2:N], nlam, bw, pos, z0, lambda.max, seg.ends, except=signal)

  R <- count_clusters_merged(x=x, z=z, z0=z0, seg.ends=seg.ends)

  fdr <- (10^(R[,1]))*p/rowSums(R[, 2:(nseg + 1), drop=FALSE])

  clust <- list()
  zsel <- matrix(ncol=nseg+1, nrow=length(level))
  for(j in 1:length(level)){
    if(any(fdr <= level[j])){
      l <- max(R[,1][fdr <= level[j]])
      ix <- which(R[,1]==l)
      zsel[j, ] <- as.numeric(z[ix,])
      clust[[j]] <- name_clusters_merged(x=x, z=rep(zsel[j,][-1], d), z0=z0)
    }else{
      clust[[j]] <-Intervals()
    }
  }
  R <- list("R" = R,  "bw"=bw, "clust"=clust, "x"=x,
            "z"=z, "z0"=z0, "zsel"=zsel)
  return(R)
}

