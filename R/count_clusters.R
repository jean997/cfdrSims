#'@import intervals
#'@import MASS

#'@title Return clusters at a given level
#'@description Find indices of clusters  given a threshold and level for merging
#'@param x Data (already smoothed)
#'@param z Threshold
#'@param z0 threshold for merging - separate regions at level
#'z that are subsets of a single region at level z0 will be merged.
#'@return Intervals object giving clusters
#'@export
name_clusters_merged <- function(x, z, z0){
  q0 <-rle( abs(x) > z0 )
  p0 <- length(q0$lengths)
  starts0 <- c(1, cumsum(q0$lengths)[-p0]+1)[q0$values]
  stops0 <- (cumsum(q0$lengths))[q0$values]
  q0I <- Intervals(cbind(starts0, stops0))

  stopifnot(length(z)==1)

  q <-rle( abs(x) > z )
  p <- length(q$lengths)
  starts <- c(1, cumsum(q$lengths)[-p]+1)[q$values]
  stops <- (cumsum(q$lengths))[q$values]
  qI <- Intervals(cbind(starts, stops))

  f <- unlist(interval_overlap(qI, q0I))
  stopifnot(length(f) == length(starts))

  t <- table(f)
  if(any(t > 1)){
    qI <- interval_merge(qI, f)
  }
  return(qI)
}

#'@title Count clusters
#'@description Count clusters  given a threshold and level for merging
#'@param x Data (already smoothed)
#'@param z Threshold
#'@param z0 threshold for merging - separate regions at level
#'z that are subsets of a single region at level z0 will be merged.
#'@return A vector the length of z giving the number of clusters at each level
#'@export
count_clusters_merged <- function(x, z, z0){
  q0 <-rle( abs(x) > z0 )
  p0 <- length(q0$lengths)
  starts0 <- c(1, cumsum(q0$lengths)[-p0]+1)[q0$values]
  stops0 <- (cumsum(q0$lengths))[q0$values]
  q0M <- cbind(starts0, stops0)
  m <- apply(q0M, MARGIN=1, FUN=function(s){ max(abs(x)[s[1]:s[2]])})
  clust_num <- unlist(lapply(z, FUN=function(t){sum(m > t)}))
  return(clust_num)
}


interval_merge <- function(obj, labs){
  n <- nrow(obj)
  t <- table(labs)
  if(all(t==1)) return(obj)
  keep.idx <- match(unique(labs), labs)
  m.names <- names(t[t > 1])
  for(m in m.names){
    strt <- min(obj[labs==m,1])
    stp <- max(obj[labs==m,2])
    i <- match(m, labs)
    obj[i,] <- c(strt, stp)
  }
  obj <- obj[keep.idx,]
  return(obj)
}
