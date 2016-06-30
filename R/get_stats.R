
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
get_stats_pois_binary <- function(dat, labs, perm_labs, s0=0){
  z1 <- pois_regression(Y=dat, labs=labs, s0=s0)
  if(is.null(perm_labs)) return(z1)
  z <- apply(perm_labs, MARGIN=2, FUN=function(l){
							 pois_regression(Y=dat,labs=l, s0=s0)
        })
  Z = cbind( z1, z)
  return(Z)
}

pois_stats2 <- function(Y, labs, s0=0){
  B <- apply(Y, MARGIN=1, FUN=function(y){
    f <- glm(y~labs, family="quasipoisson")
    b1 <- summary(f)$coefficients[2, 1]
    s <- summary(f)$coefficients[2, 2] + s0
    if(is.na(b1/s)) return(0)
    return(b1/s)
  })
  return(B)
}

get_stats_pois_continuous <- function(dat, labs, perm_labs, s0=0){
  z1 <- pois_stats2(Y=dat, labs=labs, s0=s0)
  if(is.null(perm_labs)) return(z1)
  z <- apply(perm_labs, MARGIN=2, FUN=function(l){
    pois_stats2(Y=dat,labs=l, s0=s0)
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
  z1 <- huber_stats2(Y=dat, labs=labs, s0=s0, maxit=maxit)
  if(is.null(perm_labs)) return(z1)
  z <- apply(perm_labs, MARGIN=2, FUN=function(l){
    huber_stats2(dat, labs=l, s0=s0, maxit=maxit)
  })
  Z = cbind( z1, z)
  return(Z)
}

#'Calculate Huber statsitic for simple linear regression usnig rlm
#'@description Calculate two sample Huber
#'@param Y matrix (p x n)
#'@param labs 0s and 1s
#'@param k Threshold for huber estimator in multiples of scale parameter.
#'@param s0 Additional variance
#'@param maxit Maixum iterations to pass to rlm.
#'@return Vector of test statistics.
#'@export
huber_stats2 <- function(Y, labs, k=1.345, s0=0, maxit=20){
  B <- apply(Y, MARGIN=1, FUN=function(y){
    f <- rlm(y~labs, psi=psi.huber, k=k, scale.est="Huber", maxit=maxit)
    b1 <- summary(f)$coefficients[2, 1]
    s <- summary(f)$coefficients[2, 2] + s0
    if(is.na(b1/s)) return(0)
    return(b1/s)
  })
  return(B)
}

#'Just like huber_stats2 but returns coefficient and sd estimates
#'@description Calculate two sample Huber
#'@param Y matrix (p x n)
#'@param labs 0s and 1s
#'@param k Threshold for huber estimator in multiples of scale parameter.
#'@param s0 Additional variance
#'@param maxit Maixum iterations to pass to rlm.
#'@return 2 by p matrix. Top row is coefficient estimate. Bottom row is sd estimates.
#'@export
huber_helper <- function(Y, labs, k=1.345, maxit=20){
  B <- apply(Y, MARGIN=1, FUN=function(y){
    f <- rlm(y~labs, psi=psi.huber, k=k, scale.est="Huber", maxit=maxit)
    b1 <- summary(f)$coefficients[2, 1]
    s <- summary(f)$coefficients[2, 2]
    return(c(b1, s))
  })
  return(B)
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

get_stats_lm <- function(dat, labs, perm_labs, s0=0){
  z1 <- lm_stats(Y=dat, labs=labs, s0=s0)
  if(is.null(perm_labs)) return(z1)
  z <- apply(perm_labs, MARGIN=2, FUN=function(l){
    lm_stats(Y=dat,labs=l, s0=s0)
  })
  Z = cbind( z1, z)
  return(Z)
}


lm_stats <- function(Y, labs, s0=0){
  B <- apply(Y, MARGIN=1, FUN=function(y){
    f <- lm(y~labs)
    b1 <- summary(f)$coefficients[2, 1]
    s <- summary(f)$coefficients[2, 2] + s0
    if(is.na(b1/s)) return(0)
    return(b1/s)
  })
  return(B)
}
