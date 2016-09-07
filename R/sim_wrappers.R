#'@import dplyr


#'@export
run_bin_aoas <- function(seed, prefix, n, type.sequence,
                         n.seg=c(2, 4), auto.min.length = c(50, 100, 200),
                         sample.size=c(20, 20), random.peak.loc=TRUE,
                        n.perms=500, waveQTL_loc="~/.local/bin/WaveQTL"){
  file.start <- paste0(prefix, "_", n)
  x <- rep(c(0, 1), sample.size)
  #Null functions
  g1 = function(x){return(list("ht"=5, "assoc"=0))}
  g2 <- function(x){return(list("ht"=rexp(n = 1, rate=1/5)+1.5, "assoc"=0))}
  g3 <- function(x){return(list("ht"=sample(c(1.5, 6), size=1, prob = c(0.8, 0.2)),
                                "assoc"=0))}

  #Associated functions
  g4 <- function(x){return( list("ht"=5 + x , "assoc"=1))}
  g5 <- function(x){return(list("ht"=rexp(n = 1, rate=1/5)+1.5+(2*x), "assoc"=1))}
  g6 <- function(x){return(list("ht"=sample(c(1.5, 6), size=1, prob = c(0.9-(0.25*x), 0.1 + 0.25*x)),
                                "assoc"=1))}

  pk.ht.funcs = c(g1, g2, g3, g4, g5, g6)

  R <- cfdr_sims3(x, pk.ht.funcs, type.sequence,
                  n.seg=n.seg, auto.min.length =auto.min.length,
                  seed=seed, n.perms=n.perms, s0=rep(0.05, 3),
                  level=c(0.02, 0.05, 0.1, 0.2),
                  save.data=TRUE, huber.maxit=50,
                  file.name=paste0(file.start, "_fret.RData"),
                  random.peak.loc=random.peak.loc, min.peak.sep=65)
  run_win_tests(file.start, waveQTL_loc)
}

#'@export
run_win_tests <- function(file.start, waveQTL_loc, s0=c(0, 0, 0),
                          naive.bw=c(32, 64), informed.bw=c(32, 64)){
  R <- getobj(paste0(file.start, "_fret.RData"))
  p <- dim(R$dat)[1]

  #Naive window
  for(bw in naive.bw){
    n <- paste0("naive", bw)
    k <- floor(p/bw)
    strt <- floor((p - (k*bw))/2) + 1
    strts <- ((0:(k-1))*bw) + strt
    stps <- strts + bw -1
    wins <- cbind(strts, stps)
    w_test <- window_test(wins, dat=R$dat,
                            pos=1:p, x=R$x, signal=R$signal$signal, s0=s0)
    save(w_test, file=paste0(file.start, "_", n, "_tests.RData"))
    if(log(bw, 2)==floor(log(bw, 2))){
      #w_waveqtl <-run_waveQTL(wins, dat=R$dat, x=R$x, signal=R$signal$signal,
       #                       waveQTL_loc=waveQTL_loc)
      #save(w_waveqtl, file=paste0(file.start, "_", n, "_waveqtl.RData"))
    }
  }
  #Informed windows
  for(bw in informed.bw){
    n <- paste0("informed", bw)
    wins <- R$signal$peaks
    d <- wins[,2]-wins[,1]+1
    needed <- bw - d
    left <- pmin(wins[,1]-1,  sapply(needed/2, FUN=floor))
    right <- needed-left
    wins[,1] <- wins[,1]-left
    wins[,2] <- wins[,2] + right
    wins <- as.matrix(interval_union(Intervals(wins)))
    d <- wins[,2]-wins[,1]+1
    w_test <- window_test(wins, dat=R$dat,
                          pos=1:p, x=R$x, signal=R$signal$signal, s0=s0)
    save(w_test, file=paste0(file.start, "_", n, "_tests.RData"))
    if(all(log(d, 2)==floor(log(d,2)))){
      #w_waveqtl <-run_waveQTL(wins, dat=R$dat, x=R$x, signal=R$signal$signal,
      #                      waveQTL_loc=waveQTL_loc)
      #save(w_waveqtl, file=paste0(file.start, "_", n, "_waveqtl.RData"))
    }
  }
}

#'@export
run_deseq2 <- function(file.start, naive.bw=c(32, 64), informed.bw=c(32, 64)){
  R <- getobj(paste0(file.start, "_fret.RData"))
  p <- dim(R$dat)[1]


  for(bw in naive.bw){
    n <- paste0("naive", bw)
    k <- floor(p/bw)
    strt <- floor((p - (k*bw))/2) + 1
    strts <- ((0:(k-1))*bw) + strt
    stps <- strts + bw -1
    wins <- cbind(strts, stps)
    w_deseq2 <-deseq2_test(wins, R$dat, 1:p, R$x, R$signal$signal)
    save(w_deseq2, file=paste0(file.start, "_", n, "_deseq2.RData"))
  }
  for(bw in informed.bw){
    n <- paste0("informed", bw)
    wins <- R$signal$peaks
    d <- wins[,2]-wins[,1]+1
    needed <- bw - d
    left <- pmin(wins[,1]-1,  sapply(needed/2, FUN=floor))
    right <- needed-left
    wins[,1] <- wins[,1]-left
    wins[,2] <- wins[,2] + right
    wins <- as.matrix(interval_union(Intervals(wins)))
    w_deseq2 <-deseq2_test(wins, R$dat, 1:p, R$x, R$signal$signal)
    save(w_deseq2, file=paste0(file.start, "_", n, "_deseq2.RData"))
  }
}
