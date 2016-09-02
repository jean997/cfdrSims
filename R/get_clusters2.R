#'@import fret


#'@title Count clusters at several levels of z and set threshold
#'@description Count clusters at several levels of z and set threshold. This function is mostly useful
#'for simulations when all of the smoothed statistics are in one data frame. For larger problems see the vignette.
#'@param smoothed.stats Matrix of smoothed stats.
#'@param pos Length p vector of positions
#'@param zmin Vector of length 1 or 2
#'@param level Threshold at which to control FDR.
#'@param segment.bounds matrix with two columns giving segment bounds
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
get_clusters2 <- function(smoothed.stats, pos, zmin, z0=0.3*zmin,
                         level=c(0.02, 0.05, 0.1, 0.2), segment.bounds=NULL){

  stopifnot(length(z0)==length(zmin))
  stopifnot(length(z0) %in% c(1, 2))
  if(length(z0)==2){
    z0 <- sort(z0, decreasing = TRUE)
    zmin <- sort(zmin, decreasing=TRUE)
  }

  N <- ncol(smoothed.stats)
  p <- dim(smoothed.stats)[1]
  s <- length(zmin)
  if(s==1) signed <- FALSE
    else signed=TRUE

  if(is.null(segment.bounds)) segment.bounds <- matrix(c(min(pos), max(pos)), ncol=2)

  stopifnot(ncol(segment.bounds)==2)
  stopifnot(all(segment.bounds[,1] < segment.bounds[,2]))
  K <- nrow(segment.bounds)
  stopifnot(all(segment.bounds[-K, 2] < segment.bounds[-1, 1]))

  d <- segment.bounds[,2]-segment.bounds[,1] + 1

  max1.list <- list()
  perm.maxes.list <- list()
  for(i in 1:nrow(segment.bounds)){
    ix1 <- max(which(pos <= segment.bounds[i,1]))
    ix2 <- min(which(pos >= segment.bounds[i,2]))
    max1.list[[i]] <- mxlist(ys=smoothed.stats[ix1:ix2,1], z0=z0, zmin=zmin)
    pm <- apply(smoothed.stats[ix1:ix2, 2:N], MARGIN=2, FUN=function(ys){
      mxlist(ys, z0, zmin)
    })
    pm <- sort(unlist(pm), decreasing = TRUE)
    perm.maxes.list[[i]] <- lamtab(pm, zmin, d[i], N-1)
  }

  R <- fret_choose_z2(max1.list, perm.maxes.list, nbp=d, signed=signed)

  clust <- list()
  zsel <- array(dim=c(s, length(level), K+1))
  dimnames(zsel) = list(1:s, level, names(R$z))
  for(j in 1:length(level)){
    if(any(R$fdr <= level[j])){
      ix <- max(which(R$fdr <= level[j]   ))
      zsel[1, j, ] <- as.numeric(R$z[ix,])
      if(s==2){
        zsel[2, j, ] <- as.numeric(R$zneg[ix,])
        cpos <- name_clusters_merged(x=smoothed.stats[,1],
                                     z=rep(as.numeric(zsel[1, j,-1]), d),
                                     z0=z0[1],sgn=TRUE)
        cneg <- name_clusters_merged(x=smoothed.stats[,1],
                                     z=rep(as.numeric(zsel[2, j,-1]), d),
                                     z0=z0[2],sgn=TRUE)
        clust[[j]] <- interval_union(cpos, cneg)
      }else{
        clust[[j]] <- name_clusters_merged(x=smoothed.stats[,1], z=rep(zsel[1, j,-1], d), z0=z0)
      }
    }else{
      clust[[j]] <-Intervals()
    }
  }
  R[["clust"]] <- clust
  R[["zsel"]] <- zsel
  R[["zmin"]] <- zmin
  R[["z0"]] <- z0
  return(R)
}

