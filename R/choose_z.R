


choose_z <- function(lhat, R, level){

  nseg <- ncol(lhat)
  l <- rowMeans(lhat)
  r <- rowMeans(R)
  if(all(l/r > level)){
    idx <- rep(Inf, nseg)
    f <- lhat/R
    f[lhat==0 & R ==0] <- Inf
    if(all(f > level)){
      return(idx)
    }else{
      m <- max(R[ f < level])
      ix <- which(R==m & f < level, arr.ind = TRUE)
      j <- which.min(ix[,1])
      idx[ix[j, 2]] <- ix[j, 1]
      idx <- cbind(idx, 1:nseg)
    }

  }else{
    f <- l/r
    m <- max(r[f < level])
    ix <- which(r==m & f < level)
    idx <- cbind(rep(min(ix), nseg), 1:nseg)
  }
  if(nseg==1) return(idx)
  changed <- rep(1, nseg)
  j <- 1
  while(any(changed > 0)){
    cat(j, "..")
    l  <- sum(lhat[idx[-j,]])
    r <- sum(R[idx[-j,]])
    ll <- l + lhat[,j]
    rr <- r + R[,j]
    ff <- ll/rr
    m <- max(rr[ff < level])
    ix <- which(rr==m & ff < level)
    if(min(ix)==idx[j, 1])changed[j] <- 0
      else changed[j] <- 1
    idx[j, 1] <- min(ix)
    j <- (j+1) %% nseg
    if(j==0) j <- nseg
  }
  r <- sum(R[idx])
  if(r == 0) idx[,1] <- Inf
  return(idx)

}
