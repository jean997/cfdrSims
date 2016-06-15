

ksmooth_0 <- function(x, y, bandwidth, x.range=NULL){
  n <- length(y)

  y.out <- sapply(x, FUN=function(xx){
    sum(y[ x <= (xx+bandwidth/2) & x >= (xx - bandwidth/2)])/(bandwidth+1)
  })
  if(!is.null(x.range)) y.out <- y.out[x.range[1]:x.range[2]]
  return(y.out)
}
