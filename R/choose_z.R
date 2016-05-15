
choose_z_even <- function(perm.stats, nlam, bw, pos, z0,
                          lambda.max=NULL, seg.ends=NULL, except=NULL){
  if(is.null(seg.ends)){
    seg.ends <- c(length(x))
  }
  nseg <- length(seg.ends)
  nperm <- ncol(perm.stats)
  perm.smooth <- apply(perm.stats, MARGIN=2, FUN=function(x){
    ksmooth(pos, x, bandwidth=bw, x.points=pos)$y
  })
  if(!is.null(except)){
    for(j in 1:nrow(except)) {
      s <- except[j, 1]; p <- except[j, 2]
      perm.smooth[s:p, ] <- 0
    }
  }

  mx <- list()
  max.lam <- Inf; min.lam <- -Inf; strt <- 1
  for(i in 1:nseg){
    cat(i, " ")
    maxes <- apply(perm.smooth[strt:seg.ends[i],], MARGIN=2, FUN=function(xs){
      q0 <-rle( abs(xs) > z0 )
      p0 <- length(q0$lengths)
      starts0 <- c(1, cumsum(q0$lengths)[-p0]+1)[q0$values]
      stops0 <- (cumsum(q0$lengths))[q0$values]
      sapply(1:length(starts0), FUN=function(j){ max(abs(xs)[starts0[j]:stops0[j]])})
    })
    m <- sort(unlist(maxes), decreasing=TRUE)
    mx[[i]] <- cbind(m, (1:length(m))/(nperm*(seg.ends[i]-strt + 1)))
    max.lam <- min(max.lam, max(log10(mx[[i]][,2])))
    min.lam <- max(min.lam, min(log10(mx[[i]][,2])))
    strt <- seg.ends[i] + 1
  }
  cat("\n")
  if(!is.null(lambda.max)) max.lam <- log10(lambda.max)

  lams <- seq(min.lam, max.lam, length.out=nlam)
  z <- sapply(mx, FUN=function(m){
   approx(y=m[,1], x=log10(m[,2]), xout=lams)$y
  })
  z <- data.frame(cbind(lams, z))
  names(z) <- c("lambda", paste0("z", 1:nseg))
  return(z)
}



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
    f[l==0 & r==0] <- Inf
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
    }else if(sum(ff <= level)==0){ #This only happens for numerical reasons
      changed[j] <- 0
    }else{
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
