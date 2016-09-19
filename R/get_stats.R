
#'Calculate two sample Quasi-Poisson regression statistic (binary predictor)
#'@description Calculate two sample Quasi-Poisson regression statsitic
#'@param Y matrix of data (p x n)
#'@param x trait value (0s and 1s)
#'@param s0 Additional variance
#'@param zero.val Value to use if total in one group is zero
#'@return 3 by p matrix
#'@export
qpois_stats_binary <- function(Y, x, s0, zero.val=1e-11){
  B <- apply(Y, MARGIN=1, FUN=function(y){
    c0 <- max(zero.val, sum(y[x==0]))
    c1 <- max(zero.val, sum(y[x==1]))
    n0 <- sum(x==0)
    n1 <- sum(x==1)
    b1 <- log(c1/n1)-log(c0/n0)
    mu = rep(c0/n0, n0+n1)
    mu[x==1] = c1/n1
    phi = 1/(n0 + n1 - 2)* sum( (y-mu)^2/mu)
    s = sqrt(phi)*sqrt(1/c0 + 1/c1)
    return(c(b1, s, b1/(s+s0)))
  })
  return(B)
}


#'Calculate two sample Poisson regression statsitic for any type of predictor
#'@description Calculate Quasi-Poisson regression statsitic
#'@param Y matrix of data (p x n)
#'@param x trait value (0s and 1s)
#'@param s0 Additional variance
#'@return 3 by p matrix
#'@export
qpois_stats_continuous <- function(Y, x, s0){
  B <- apply(Y, MARGIN=1, FUN=function(y){
    f <- glm(y~x, family="quasipoisson")
    b1 <- summary(f)$coefficients[2, 1]
    s <- summary(f)$coefficients[2, 2]
    return(c(b1, s, b1/(s+s0)))
  })
  return(B)
}


#'Calculate two sample t-statsitic (binary predictor)
#'@description Calculate two sample t-test statsitic
#'@param Y matrix of data (p x n)
#'@param x trait value (0s and 1s)
#'@param s0 Additional variance
#'@return 3 by p matrix
#'@export
t_stats <- function(Y, x, s0){
  n0 <- sum(x==0); n1 <- sum(x==1)
  B <- apply(Y, MARGIN=1, FUN=function(y){
    m0 <- sum(y[x==0])/n0
    v0 <- var(y[x==0])
    m1 <- sum(y[x==1])/n1
    v1 <- var(y[x==1])
    b1 <- m1 -m0
    s <- sqrt(v0/n0 + v1/n1)
    return(c(b1, s, b1/(s+s0)))
  })
  return(B)
}


#'Calculate two sample linear regression test statistics
#'@description Calculate linear regression test statsitic
#'@param Y matrix of data (p x n)
#'@param x trait value (0s and 1s)
#'@param s0 Additional variance
#'@return 3 by p matrix
#'@export
lm_stats <- function(Y, x, s0){
  B <- apply(Y, MARGIN=1, FUN=function(y){
    f <- lm(y~x)
    b1 <- summary(f)$coefficients[2, 1]
    s <- summary(f)$coefficients[2, 2]
    return(c(b1, s, b1/(s+s0)))
  })
  return(B)
}
