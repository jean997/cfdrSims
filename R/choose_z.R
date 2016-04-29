


choose_z <- function(lhat, R, level){

  nseg <- ncol(lhat)
  l <- rowMeans(lhat)
  r <- rowMeans(R)
  if(all(l/r > level)){
    idx <- cbind(rep(Inf, nseg), 1:nseg)
    f <- lhat/R
    f[lhat==0 & R ==0] <- Inf
    if(all(f > level)) return(idx)
  }else{
    f <- l/r
    m <- max(r[f <= level])
    ix <- which(r==m & f <= level)
    idx <- cbind(rep(min(ix), nseg), 1:nseg)
  }
  if(nseg==1) return(idx)
  changed <- rep(1, nseg)
  j <- 1
  while(any(changed > 0)){
    cat(j, "..")
    inf_ix <- which(!is.finite(idx[,1]))
    l  <- sum(lhat[idx[-c(inf_ix, j),, drop=FALSE]])
    r <- sum(R[idx[-c(inf_ix, j),, drop=FALSE]])
    ll <- l + lhat[,j]
    rr <- r + R[,j]
    ff <- ll/rr
    ff[ll==0 & rr==0] <- Inf
    if(all(ff > level) & idx[j, 1]==Inf){
      changed[j] <- 0
    }else{
      if(sum(ff <= level)==0) cat(j, "!!!\n")
      m <- max(rr[ff <= level])
      ix <- which(rr==m & ff <= level)
      if(min(ix)==idx[j, 1])changed[j] <- 0
        else changed[j] <- 1
      idx[j, 1] <- min(ix)
    }
    j <- (j+1) %% nseg
    if(j==0) j <- nseg
  }
  inf_ix <- which(!is.finite(idx[,1]))
  if(length(inf_ix)==0) r <- sum(R[idx])
    else r <- sum(R[idx[-c(inf_ix),, drop=FALSE] ])
  if(r == 0) idx[,1] <- Inf
  return(idx)

}
