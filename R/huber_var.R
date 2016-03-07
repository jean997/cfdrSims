#'Calculate variance of the huber estimator
#'@description Calculate variance of the huber estimator using influence functions
#'@param x Data
#'@param muH Estimate of true huber mean
#'@param k Absolute threshold
#'@return variance
#'@export
huber_var <- function(x, muH, k){
  n <- length(x)
  ifs <- x-muH
  ifs[(x-muH) < -k] <- -k
  ifs[(x-muH) > k] <- k
  C <- sum(abs(x-muH)<=k)/n
  ifs <- ifs/C
  return(sum(ifs^2)/n^2)
}


#'Calculate two sample Huber statsitic
#'@description Calculate two sample Huber
#'@param Y matrix (p x n)
#'@param labs 0s and 1s
#'@param k Threshold for huber estimator in multiples of scale parameter.
#'@param scale Vector of length p giving scale parameter.
#'Will be calculated if not provided.
#'@param s0 Additional variance
#'@param zero.val Value to use if total in one group is zero
#'@return Vector of statistics length p
#'@export
huber_stats <- function(Y, labs, k=1.345, scale=NULL, s0=0){
  if(is.null(scale)){
    scale <- huber_scale_prop2(Y, labs, k)
  }
  if(length(scale)==1) scale <- rep(scale, nrow(Y))
  Yt <- cbind(scale, Y)
  B <- apply(Yt, MARGIN=1, FUN=function(y){
    s <- y[1]; y <- y[-1]
    b0 <- hubers(y[labs==0], k=k*s, s=1)$mu
    b1 <- hubers(y[labs==1], k=k*s, s=1)$mu
    v0 <- huber_var(y[labs==0], b0, k*s)
    v1 <- huber_var(y[labs==1], b1, k*s)
    v <- v0 + v1
    return((b1-b0)/(sqrt(v)+s0))
  })
  return(B)
}

huber_scale_mad <- function(Y, labs, k2=1.345){
  R <- matrix(nrow=nrow(Y), ncol=ncol(Y))
  m1 <- apply(Y[,labs==0], MARGIN=1, FUN=median)
  m2 <- apply(Y[,labs==1], MARGIN=1, FUN=median)
  R[,labs==0] <- Y[,labs==0]-m1
  R[,labs==1] <- Y[,labs==1]-m2

  s1 <- apply(R, MARGIN=1, FUN=function(resids){mad(resids, 0)})
  s1[s1==0] <- min(s1[s1>0])/2
  return(s1)
}

huber_scale_prop2 <- function(Y, labs, k2=1.345){

  R <- matrix(nrow=nrow(Y), ncol=ncol(Y))
  m1 <- apply(Y[,labs==0], MARGIN=1, FUN=median)
  m2 <- apply(Y[,labs==1], MARGIN=1, FUN=median)
  R[,labs==0] <- Y[,labs==0]-m1
  R[,labs==1] <- Y[,labs==1]-m2

  s1 <- apply(R, MARGIN=1, FUN=function(resids){mad(resids, 0)})
  s1[s1==0] <- min(s1[s1>0])/2
  n1 <- ncol(R)-2
  theta <- 2*pnorm(k2) - 1
  gamma <- theta + k2^2 * (1 - theta) - 2 * k2 * dnorm(k2)
  s <- apply(cbind(s1, R), MARGIN=1, FUN=function(resids){
    scale <- resids[1]
    resids <- resids[-1]
    sqrt(sum(pmin(resids^2, (k2 * scale)^2))/(n1*gamma))
  })
  return(s)
}

huber_stats2 <- function(Y, labs, k=1.345, s0=0){
  B <- apply(Y, MARGIN=1, FUN=function(y){
    f <- rlm(y~labs, psi=psi.huber, k=k, scale.est="Huber")
    b1 <- summary(f)$coefficients[2, 1]
    s <- summary(f)$coefficients[2, 2] + s0
    return(b1/s)
  })
  return(B)
}
