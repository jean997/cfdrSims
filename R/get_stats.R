#set.seed(17897563)
#perms <- replicate(n=500, expr = {
#   sample( rep(c(0, 1), each=10), size=20, replace=FALSE)
#  })
#save(perms, file="perm_labs.RData")

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
get_stats_pois <- function(dat, labs, perm_labs, s0=0.1){
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
get_stats_huber<- function(dat, labs, perm_labs, s0=0.1){
  z1 <- huber_stats(Y=dat, labs=labs, s0=s0)
  z <- apply(perm_labs, MARGIN=2, FUN=function(l){
    huber_stats(dat, labs=l, s0=s0)
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
  z <- apply(perm_labs, MARGIN=2, FUN=function(l){
    t_stats(dat, labs=l, s0=s0)
  })
  Z = cbind(z1, z)
  return(Z)
}

t_stats <- function(dat, labs, s0){
  stats <- apply(dat, MARGIN=1, FUN=function(x){
    m <- as.vector(by(data=x, INDICES = labs, FUN=mean))
    v <- as.vector(by(data=x, INDICES = labs, FUN=var))
    (m[1]-m[2])/(sqrt(v[1]/sum(labs==0) + v[2]/sum(labs==1)) + s0)
  })
  return(stats)
}
