
get_signal <- function(type.sequence, width=200, peak.base=20, bandwidth=20){
  r <- length(type.sequence)
  q <- sum(type.sequence %in% 4:6)
  p <- width * length(type.sequence)

  midpoint <- width/2
  rad <- peak.base/2 + bandwidth/2
  rng <- midpoint + c(-1, 1)*rad

  rI <- Intervals(cbind(width*((1:r)-1) + 1, width*(1:r)))
  if(q==0){
    signal <- Intervals()
    R <- list("signal"=signal, "regions"=rI)
    return(R)
  }

  signal <- matrix(nrow=q, ncol=2)
  i <- 1; ct <- 1
  for(i in 1:r){
    if(type.sequence[i] %in% 4:6){
      signal[ct,] <- rng + (i-1)*width
      ct <- ct + 1
    }
  }
  signal <- Intervals(signal)
  R <- list("signal"=signal, "regions"=rI)
  return(R)
}
