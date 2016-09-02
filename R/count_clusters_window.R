

count_clusters_window <- function(x, z, z0, pos, win.size){
  stopifnot(length(x)==length(pos))
  stopifnot(all(is.integer(pos)))
  clust  = name_clusters_merged(x, z, z0)
  if(nrow(clust)==0) return(rep(0, length(x)))
  clust[,1] <- pos[clust[,1]]
  clust[,2] <- pos[clust[,2]]

  cnum <- rep(0, length(x))
  for(i in 1:nrow(clust)){
    ix <- which(pos >= clust[i, 1]-win.size & pos <= clust[i, 2]+ win.size)
    cnum[ix] <- cnum[ix] + 1
  }
  return(cnum)
}
