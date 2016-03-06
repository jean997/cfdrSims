#' Sample data
#' @description Sample data given a sample size and a vector of types
#' @param type.sequence Vector of integers in 1:6 giving the type of each subregion.
#' @param sample.size Vector of length 2 giving sample sizes
#' @return A list with elements:
#' \describe{
#' \item{\code{dat}}{ Matrix of observed data (200*w) x sum(sample.size) where w is length(type.sequence)}
#' \item{\code{means}}{ Matrix of sample specific profiles.}
#' }
#' @export
sample_data <- function(type.sequence, sample.size, type.def=NULL){
  stopifnot(length(sample.size)==2)

  P <- define_profiles()
  if(is.null(type.def)) type.def <- define_types()
  z <- nrow(P)
  w <- length(type.sequence)
  dat <- means <- matrix(nrow=w*z, ncol=sum(sample.size))

  for(i in 1:sample.size[1]){
    idx <- 1
    for(j in 1:w){
      p <- P[,sample(1:4, size=1, prob=type.def$p1[type.sequence[j],])]
      dat[idx:(idx+z-1),i] <- rpois(n=z, lambda = p)
      means[idx:(idx+z-1), i] <- p
      idx <- idx + z
    }
  }
  for(i in 1:sample.size[2]){
    idx <- 1
    for(j in 1:w){
      p <- P[,sample(1:4, size=1, prob=type.def$p2[type.sequence[j],])]
      dat[idx:(idx+z-1),i+sample.size[1]] <- rpois(n=z, lambda = p)
      means[idx:(idx+z-1), i+sample.size[1]] <- p
      idx <- idx + z
    }
  }
  return(list("dat"=dat, "means"=means))
}
