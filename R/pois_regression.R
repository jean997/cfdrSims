#'Calculate two sample Poisson regression statsitic
#'@description Calculate two sample Poisson regression statsitic
#'@param Y matrix of data (p x n)
#'@param labs 0s and 1s
#'@param s0 Additional variance
#'@param zero.val Value to use if total in one group is zero
#'@return A length p vector of statistics
#'@export
pois_regression <- function(Y, labs, s0, zero.val=1e-11){
  B <- apply(Y, MARGIN=1, FUN=function(y){
    c0 <- max(zero.val, sum(y[labs==0]))
    c1 <- max(zero.val, sum(y[labs==1]))
    n0 <- sum(labs==0)
    n1 <- sum(labs==1)
    beta1 <- log(c1/n1)-log(c0/n0)
    mu = rep(c0/n0, n0+n1)
    mu[labs==1] = c1/n1
    phi = 1/(n0 + n1 - 2)* sum( (y-mu)^2/mu)
    s1 = sqrt(phi)*sqrt(1/c0 + 1/c1)
    return(beta1/(s1+s0))
  })
  return(B)
}

