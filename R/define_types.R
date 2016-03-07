define_types <- function(){
  p1 <- p2 <- matrix(nrow=6, ncol=4)

  p1[1,] <- p2[1,] <- c(0, 0, 1, 0)
  p1[2,] <- p2[2,] <- c(0.2, 0.2, 0.5, 0.1)
  p1[3,] <- p2[3,] <- c(0.7, 0.1, 0.25, 0.1)

  p1[4,] <- c(1, 0, 0, 0)
  p2[4,] <- c(0.6, 0.4, 0, 0)

  p1[5,] <- c(0.2, 0.2, 0.5, 0.1)
  p2[5,] <- c(0.55, 0.3, 0.1, 0.05)

  p1[6,] <- c(0.5, 0, 0.4, 0.1)
  p2[6,] <- c(0.95, 0, 0, 0.05)
  return(list("p1"=p1, "p2"=p2))
}

signal_to_noise <- function(p1, p2, h){
  m1 <- sum(p1*h)
  v1 <- sum(p1*(h^2 + h)) - m1^2
  m2 <- sum(p2*h)
  v2 <- sum(p2*(h^2 + h)) - m2^2
  return(list("mu"=m1-m2, "sigma"=sqrt(v1+v2), "snr"=abs(m1-m2)/sqrt(v1+v2)))

}
