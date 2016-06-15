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

huber_var_eq6.5 <- function(f){
  stopifnot(any(class(f)=="rlm"))
  X <- model.matrix(f)
  r <- f$residuals
  k <- f$k2*f$s
  m <- sum(-k < r & r < k)
  n <- length(r)
  p <- dim(X)[2]
  K <- 1 + (p/n)*(1-m)/m #Eq. (6.11), page 174 Huber, Robust Statistics
  psi_r <- r*pmin(1, k/abs(r))
  psi_prime_r <- as.numeric(abs(r) < k)
  return((K^2)*( (1/(n-p))* sum(psi_r^2)/(mean(psi_prime_r)^2))*solve(t(X)%*%X))
}


huber_var_eq6.6 <- function(f){
  stopifnot(any(class(f)=="rlm"))
  X <- model.matrix(f)
  r <- f$residuals
  k <- f$k2*f$s
  m <- sum(-k < r & r < k)
  n <- length(r)
  p <- dim(X)[2]
  K <- 1 + (p/n)*(1-m)/m #Eq. (6.11), page 174 Huber, Robust Statistics
  psi_r <- r*pmin(1, k/abs(r))
  psi_prime_r <- as.numeric(abs(r) < k)
  W <- matrix(nrow=p, ncol=p)
  for(j in 1:p){
    for(k in 1:p){
      W[j, k] <- 0
      for(i in 1:n){
        W[j, k] <- W[j, k] + psi_prime_r[i]*X[i, j]*X[i, k]
      }
    }
  }
  return(K*( (1/(n-p))* sum(psi_r^2)/mean(psi_prime_r))*solve(W))
}

huber_var_eq6.7 <- function(f){
  stopifnot(any(class(f)=="rlm"))
  X <- model.matrix(f)
  r <- f$residuals
  k <- f$k2*f$s
  m <- sum(-k < r & r < k)
  n <- length(r)
  p <- dim(X)[2]
  K <- 1 + (p/n)*(1-m)/m #Eq. (6.11), page 174 Huber, Robust Statistics
  psi_r <- r*pmin(1, k/abs(r))
  psi_prime_r <- as.numeric(abs(r) < k)
  W <- matrix(nrow=p, ncol=p)
  for(j in 1:p){
    for(k in 1:p){
      W[j, k] <- 0
      for(i in 1:n){
        W[j, k] <- W[j, k] + psi_prime_r[i]*X[i, j]*X[i, k]
      }
    }
  }
  return((1/K)* (1/(n-p))* sum(psi_r^2)*solve(W)%*%(t(X)%*%X)%*%solve(W))
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
  m1 <- apply(Y[,labs==0, drop=FALSE], MARGIN=1, FUN=median)
  m2 <- apply(Y[,labs==1, drop=FALSE],  MARGIN=1, FUN=median)
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

huber_helper <- function(Y, labs, k=1.345, maxit=20){
  B <- apply(Y, MARGIN=1, FUN=function(y){
    f <- rlm(y~labs, psi=psi.huber, k=k, scale.est="Huber", maxit=maxit)
    b1 <- summary(f)$coefficients[2, 1]
    s <- summary(f)$coefficients[2, 2]
    return(c(b1, s))
  })
  return(B)
}
