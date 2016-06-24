
#'Get Poisson stats
#'@description Calculate Poisson regression statistic at each position
#'@param dat Matrix of data (n x p)
#'@param labs Labels (0s and 1s)
#'@param perm_labs Permutation labels (n x n.perm)
#'@param s0 Additional variance
#'@return A matrix p x (1 + n.perm)
#'The first column is for te original labels.
#'The next columns are for the permutation labels
#'@export
get_stats_pois <- function(dat, labs, perm_labs, s0=0){
  z1 <- pois_regression(Y=dat, labs=labs, s0=s0)
  z <- apply(perm_labs, MARGIN=2, FUN=function(l){
							 pois_regression(Y=dat,labs=l, s0=s0)
        })
  Z = cbind( z1, z)
  return(Z)
}

#'Get Huber stats
#'@description Calculate Huber statistic at each position
#'@param dat Matrix of data (n x p)
#'@param labs Labels (0s and 1s)
#'@param perm_labs Permutation labels (n x n.perm)
#'@param s0 Additional variance
#'@return A matrix p x (1 + n.perm)
#'The first column is for te original labels.
#'The next columns are for the permutation labels
#'@export
get_stats_huber<- function(dat, labs, perm_labs,
                           s0=0, k=1.345, scale=NULL){
  if(is.null(scale)){
    scale <- huber_scale_prop2(dat, labs, k)
  }
  z1 <- huber_stats(Y=dat, labs=labs, s0=s0, k=k, scale=scale)
  z <- apply(perm_labs, MARGIN=2, FUN=function(l){
    huber_stats(dat, labs=l, s0=s0, k=k, scale=scale)
  })
  Z = cbind( z1, z)
  return(Z)
}

#'Get Huber stats
#'@description Calculate Huber statistic at each position
#'@param dat Matrix of data (n x p)
#'@param labs Labels (0s and 1s)
#'@param perm_labs Permutation labels (n x n.perm)
#'@param s0 Additional variance
#'@param maxit Maximum iterations for rlm
#'@return A matrix p x (1 + n.perm)
#'The first column is for te original labels.
#'The next columns are for the permutation labels
#'@export
get_stats_huber2<- function(dat, labs, perm_labs, s0=0, maxit=50){
  z1 <- cfdrSims:::huber_stats2(Y=dat, labs=labs, s0=s0, maxit=maxit)
  if(is.null(perm_labs)) return(z1)
  z <- apply(perm_labs, MARGIN=2, FUN=function(l){
    cfdrSims:::huber_stats2(dat, labs=l, s0=s0, maxit=maxit)
  })
  Z = cbind( z1, z)
  return(Z)
}

#'Get t-stats
#'@description Calculate t-statistic at each position
#'@param dat Matrix of data (n x p)
#'@param labs Labels (0s and 1s)
#'@param perm_labs Permutation labels (n x n.perm)
#'@param s0 Additional variance
#'@return A matrix p x (1 + n.perm)
#'The first column is for te original labels.
#'The next columns are for the permutation labels
#'@export
get_stats_ttest<- function(dat, labs, perm_labs, s0=0.1){
  z1 <- t_stats(dat, labs=labs, s0=s0)
  if(is.null(perm_labs)) return(z1)
  z <- apply(perm_labs, MARGIN=2, FUN=function(l){
    t_stats(dat, labs=l, s0=s0)
  })
  Z = cbind(z1, z)
  return(Z)
}

t_stats <- function(dat, labs, s0){
  n1 <- sum(labs==0); n2 <- sum(labs==1)
  m1 <- apply(dat[, labs==0], MARGIN=1, FUN=sum)/n1
  v1 <- apply(dat[, labs==0], MARGIN=1, FUN=var)
  m2 <- apply(dat[, labs==1], MARGIN=1, FUN=sum)/n2
  v2 <- apply(dat[, labs==1], MARGIN=1, FUN=var)
  return((m1-m2)/(sqrt(v1/n1 + v2/n2) + s0))
}
