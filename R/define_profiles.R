
lump <- function(strt, stp, ht, mesa.width=7, bw=5, bg.ht=1.5){
  total.width = stp-strt + 1
  w <- (total.width-mesa.width)/2
  y <- c(bg.ht + (1:w)*(ht-bg.ht)/w, rep(ht, mesa.width), ht- (1:w)*(ht-bg.ht)/w)
  y <- c(rep(bg.ht, bw+1), y, rep(bg.ht, bw+1))
  ys <- ksmooth(x=(strt-bw-1):(stp+bw+1),y= y, bandwidth = bw, x.points=strt:stp)$y
  return(ys)
}

gen_profile2 <- function(peak.starts, peak.hts, total.length,
                         peak.base=20, mesa.width=7, bw=5, bg.ht=1.5){
  y <- rep(bg.ht, total.length)
  stopifnot(length(peak.starts)==length(peak.hts))
  d <- diff(peak.starts)
  stopifnot(all(d >= peak.base))
  for(i in 1:length(peak.hts)){
    strt <- peak.starts[i]
    stp <- peak.starts[i] + peak.base -1
    y[strt:stp] <- lump(strt, stp, peak.hts[i], mesa.width, bw, bg.ht)
  }
  return(y)

}

get_signal3 <- function(pk.ht.funcs, type.sequence, peak.starts, peak.base=20){
  signal <- peaks <- cbind(peak.starts, peak.starts + peak.base -1)
  is.signal <- sapply(type.sequence, FUN=function(t){
    pk.ht.funcs[[t]](0)$assoc
  })
  signal <- signal[is.signal==1,]
  return(list("signal"=signal, "peaks"=peaks))
}
