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
  if(is.null(scale)) scale <- apply(Y, MARGIN=1, FUN=sd)
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

huber_stats2 <- function(Y, labs, k=1.345, s0=0){
  B <- apply(Y, MARGIN=1, FUN=function(y){
    f <- rlm(y~labs, psi=psi.huber, k=k, scale.est="Huber")
    b1 <- summary(f)$coefficients[2, 1]
    s <- summary(f)$coefficients[2, 2] + s0
    return(b1/s)
  })
  return(B)
}
