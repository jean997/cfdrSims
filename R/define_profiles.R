define_profiles <- function(){
  bg <- 1.5
  width <- 200
  peak.base <- 20
  hts <- c(bg, 3, 10, 60)
  profiles <- matrix(nrow=width, ncol=length(hts))
  for(i in 1:length(hts)){
    w <- peak.base/2
    w1 <- width/2 - peak.base/2
    w2 <- w1 + peak.base/2
    w3 <- w2 + peak.base/2
    w4 <- w3 + w1
    profiles[1:w1,i] <- bg
    profiles[(w3+1):w4, i] <- bg
    profiles[(w1+1):w2, i] <- bg + (1:w)*(hts[i]-bg)/w
    profiles[(w2+1):w3, i] <- hts[i]- (1:w)*(hts[i]-bg)/w
  }
  return(profiles)

}

lump <- function(strt, stp, ht, mesa.width, bw=0, bg.ht=1.5){
  total.width = stp-strt + 1
  w <- (total.width-mesa.width)/2
  y <- c(bg.ht + (1:w)*(ht-bg.ht)/w, rep(ht, mesa.width), ht- (1:w)*(ht-bg.ht)/w)
  y <- c(rep(bg.ht, bw+1), y, rep(bg.ht, bw+1))
  ys <- ksmooth(x=(strt-bw-1):(stp+bw+1),y= y, bandwidth = bw, x.points=strt:stp)$y
  return(ys)
}

gen_profile <- function(ht, bg.ht=1.5, total.width=200, peak.base=20, mesa.width=7, bw=5, n.sep=2){
  peak.total.base <- peak.base*length(ht) + n.sep*(length(ht)-1)
  w1 <- floor((total.width-peak.total.base)/2)
  w2 <- total.width-peak.total.base-w1
  y <- rep(bg.ht, w1-n.sep)
  for(i in 1:length(ht)){
    y <- c(y, rep(bg.ht, n.sep))
    y <- c(y, lump(1, peak.base, ht[i], mesa.width, bw, bg.ht))
  }
  y <- c(y, rep(bg.ht, w2))
}

lump_intervals <- function(n.lumps=1, total.width=200, peak.base=20, n.sep=2){
  peak.total.base <- peak.base*n.lumps + n.sep*(n.lumps-1)
  w1 <- floor((total.width-peak.total.base)/2)
  w2 <- total.width-peak.total.base-w1
  starts <- w1 + (1:n.lumps-1)*(peak.base+n.sep)
  stops <- starts + peak.base
  return(cbind(starts, stops))
}

get_signal2 <- function(pk.ht.funcs, type.sequence, total.width=200, peak.base=20, n.sep=2){
  signal.ivls <- list()
  peak.ivls <- list()
  i <- 1
  for(f in pk.ht.funcs){
    assoc <- f(0)$assoc
    peak.ivls[[i]] <- iv <-  lump_intervals(n.lumps=length(assoc),
                                     total.width = total.width, peak.base=peak.base, n.sep=n.sep)
    if(all(assoc==0)){
      signal.ivls[[i]] <- matrix(nrow=0, ncol=2)
      next
    }
    iv <- iv[assoc==1, ]
    signal.ivls[[i]] <- iv
    i <- i+1
  }
  ct = 0
  signal <- peaks <-  matrix(nrow=0, ncol=2)
  for(t in type.sequence){
    peaks <- rbind(peaks, ct + peak.ivls[[t]])
    if(nrow(signal.ivls[[t]])==0){
      ct = ct + total.width
      next
    }
    signal <- rbind(signal, ct + signal.ivls[[t]])
    ct = ct + total.width
  }
  return(list("signal"=signal, "peaks"=peaks))
}
