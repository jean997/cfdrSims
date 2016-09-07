#'@import dplyr


#'@export
run_bin_aoas <- function(seed, prefix, n, type.sequence,
                         n.seg=c(2, 4), auto.min.length = c(50, 100, 200),
                         sample.size=c(10, 20),
                    n.perms=500, waveQTL_loc="~/.local/bin/WaveQTL"){
  file.start <- paste0(prefix, "_", n)
  x <- rep(c(0, 1), sample.size)
  g1 = function(x){return(list("ht"=5, "assoc"=0))}
  g2 <- function(x){return(list("ht"=rexp(n = 1, rate=1/5)+1.5, "assoc"=0))}
  g3 <- function(x){
    if(x==0) return( list("ht"=4, "assoc"=1))
    return(list("ht"=5, "assoc"=1))
  }
  g4 <- function(x){return(list("ht"=rexp(n = 1, rate=1/5)+1.5+x, "assoc"=1))}
  pk.ht.funcs = c(g1, g2, g3, g4)

  R <- cfdr_sims3(x, pk.ht.funcs, type.sequence,
                  n.seg=n.seg, auto.min.length =auto.min.length,
                  seed=seed, n.perms=n.perms, s0=rep(0.05, 3),
                  level=c(0.02, 0.05, 0.1, 0.2),
                  save.data=TRUE, huber.maxit=50,
                  file.name=paste0(file.start, "_fret.RData"))
  run_win_tests(file.start, p = 200*length(type.sequence), waveQTL_loc )
}







#'@export
run_bin <- function(seed, prefix, n, type.sequence, n.seg=c(2, 6), sample.size=c(15, 15),
                    n.perm=500, waveQTL_loc="~/.local/bin/WaveQTL"){
  file.start <- paste0(prefix, "_", n)
  x <- rep(c(0, 1), sample.size)
  g1 = function(x){return(list("ht"=5, "assoc"=0))}
  g2 <- function(x){return(list("ht"=rexp(n = 1, rate=1/5)+1.5, "assoc"=0))}
  g3 <- function(x){
    if(x==0) return( list("ht"=sample(c(1.5, 3, 7, 10), prob=c(0.55, 0.3, 0.1, 0.05), size=1), "assoc"=1))
    return(list("ht"=sample(c(1.5, 3, 7, 10), prob=c(0.2, 0.45, 0.25, 0.1), size=1), "assoc"=1))
  }
  g4 <- function(x){
    p1 = sample(c(1.5, 3, 7, 10), prob=c(0.55, 0.3, 0.1, 0.05), size=1)
    p2 = sample(c(1.5, 3, 7, 10), prob=c(0.2, 0.45, 0.25, 0.1), size=1)
    if(x==0) return( list("ht"=c(p1, p2), "assoc"=c(1, 1)))
    return(list("ht"=c(p2, p1), "assoc"=c(1, 1)))
  }
  g5 <- function(x){
    if(x==0) return( list("ht"=4, "assoc"=1))
    return(list("ht"=5, "assoc"=1))
  }
  g6 <- function(x){return(list("ht"=rexp(n = 1, rate=1/5)+1.5+x, "assoc"=1))}
  pk.ht.funcs = c(g1, g2, g3, g4, g5, g6)

  R <- cfdr_sims2(x, pk.ht.funcs, type.sequence, n.seg=n.seg,
                         seed=seed, n.perms=n.perm, s0=rep(0.05, 3),
                         level=c(0.02, 0.05, 0.1, 0.2),
                         save.data=TRUE, huber.maxit=50,
                         file.name=paste0(file.start, "_fret.RData"))
  run_win_tests(file.start, p = 200*length(type.sequence), waveQTL_loc )
}


#'@export
run_cont <- function(seed, prefix, n, type.sequence, n.seg=c(2, 6), sample.size=30,
                       n.perm=500, waveQTL_loc="~/.local/bin/WaveQTL"){
  file.start <- paste0(prefix, "_", n)
  x <- rbeta(n=sample.size, shape1=2, shape2=2)

  #Height Functions
  g1 = function(x){return(list("ht"=5, "assoc"=0))}
  g2 <- function(x){return(list("ht"=rexp(n = 1, rate=1/5)+1.5, "assoc"=0))}
  g3 <- function(x){
    return(list("ht"=3+4*x, "assoc"=1))
  }
  g4 <- function(x){
    p1 <- 3 + 4*x
    p2 <- sample(c(0, 3), prob=c(0.9, 0.1), size=1)
    return(list("ht"=p1+p2, "assoc"=1))
  }
  g5 <- function(x){
    p1 <- sample(c(1.5, 7), prob=c(1-x, x), size=1)
    p2 <- sample(c(0, 3), prob=c(0.9, 0.1), size=1)
    return(list("ht"=p1+p2, "assoc"=1))
  }
  g6 <- function(x){
    p1 <- c(3+4*x, 7-4*x)
    p2 <- sample(c(0, 3), prob=c(0.9, 0.1), size=2, replace=TRUE)
    return(list("ht"=p1+p2, "assoc"=c(1, 1)))
  }
  g7 <- function(x){
    p1 <-  4*x +  rexp(n=1, rate=1/4)
    return(list("ht"=p1, "assoc"=1))
  }
  pk.ht.funcs = c(g1, g2, g3, g4, g5, g6, g7)


  R <- cfdr_sims2(x, pk.ht.funcs, type.sequence, n.seg=n.seg,
                  seed=seed, n.perms=n.perm, s0=rep(0.05, 3),
                  level=c(0.02, 0.05, 0.1, 0.2),
                  save.data=TRUE, huber.maxit=50,
                  file.name=paste0(file.start, "_fret.RData"))
  R <- getobj(paste0(file.start, "_fret.RData"))
  run_win_tests(file.start, p = 200*length(type.sequence), waveQTL_loc)
}


run_win_tests <- function(file.start, p, waveQTL_loc, s0=c(0, 0, 0)){
  R <- getobj(paste0(file.start, "_fret.RData"))

  #Even Windows
  w_50e <- cbind(seq(26, p-27, by=50), seq(75, p, by=50))
  w_50e_test <- window_test(w_50e, dat=R$dat, pos=1:p,
                            x=R$x, signal=R$signal$signal, s0=s0)
  save(w_50e_test, file=paste0(file.start, "_w50e.RData"))

  #Windows around peaks
  k <- p/200
  w_50b <- cbind( 76 + 200*(1:k -1), 125 + 200*(1:k -1))
  w_50b_test <- window_test(w_50b, dat=R$dat,pos=1:p,
                            x=R$x, signal=R$signal$signal, s0=s0)
  save(w_50b_test, file=paste0(file.start, "_w50b.RData"))

  #Width 64 bins for waveQTL
  #Even

  w_64e <- cbind(seq(1, p - (p%%64), by=64), seq(64, p, by=64))
  w_64e_test <- window_test(w_64e, dat=R$dat,
                            pos=1:p, x=R$x, signal=R$signal$signal, s0=s0)
  save(w_64e_test, file=paste0(file.start, "_w64e.RData"))
  w_64e_waveqtl <-run_waveQTL(w_64e, dat=R$dat, x=R$x, signal=R$signal$signal,
                      waveQTL_loc=waveQTL_loc)
  save(w_64e_waveqtl, file=paste0(file.start, "_w64e_wave.RData"))
  #At peaks
  w_64b <- cbind( 69 + 200*(1:k -1), 132 + 200*(1:k -1))
  w_64b_test <- window_test(w_64b, dat=R$dat,
                            pos=1:p, x=R$x, signal=R$signal$signal, s0=s0)
  save(w_64b_test, file=paste0(file.start, "_w64b.RData"))
  w_64b_waveqtl <-run_waveQTL(w_64b, dat=R$dat, x=R$x, signal=R$signal$signal,
                              waveQTL_loc=waveQTL_loc)
  save(w_64b_waveqtl, file=paste0(file.start, "_w64b_wave.RData"))
}

run_deseq2 <- function(file.start){
  R <- getobj(paste0(file.start, "_fret.RData"))
  p <- dim(R$dat)[1]
  k <- p/200
  #w_64e <- cbind(seq(1, p - (p%%64), by=64), seq(64, p, by=64))
  #w_64b <- cbind( 69 + 200*(1:k -1), 132 + 200*(1:k -1))
  #w_64e_deseq2 <-deseq2_test(w_64e, R$dat, 1:p, R$x, R$signal$signal)
  #save(w_64e_deseq2, file=paste0(file.start, "_w64e_deseq2.RData"))
  #w_64b_deseq2 <-deseq2_test(w_64b, R$dat, 1:p, R$x, R$signal$signal)
  #save(w_64b_deseq2, file=paste0(file.start, "_w64b_deseq2.RData"))

  n <- "w64e"
  w <-  cbind(seq(1, p - (p%%64), by=64), seq(64, p, by=64))
  w_deseq2 <-deseq2_test(w, R$dat, 1:p, R$x, R$signal$signal)
  save(w_deseq2, file=paste0(file.start, "_", n, "_deseq2.RData"))

  n <- "w64b"
  w <-  cbind( 69 + 200*(1:k -1), 132 + 200*(1:k -1))
  w_deseq2 <-deseq2_test(w, R$dat, 1:p, R$x, R$signal$signal)
  save(w_deseq2, file=paste0(file.start, "_", n, "_deseq2.RData"))


  n <- "w32e"
  w <-  cbind(seq(1, p - (p%%32), by=32), seq(32, p, by=32))
  w_deseq2 <-deseq2_test(w, R$dat, 1:p, R$x, R$signal$signal)
  save(w_deseq2, file=paste0(file.start, "_", n, "_deseq2.RData"))

  n <- "w32b"
  w <-  cbind( 85 + 200*(1:k -1), 116 + 200*(1:k -1))
  w_deseq2 <-deseq2_test(w, R$dat, 1:p, R$x, R$signal$signal)
  save(w_deseq2, file=paste0(file.start, "_", n, "_deseq2.RData"))


  #n <- "w64o"
  #w <- cbind( 69 + 200*(1:k -1), 132 + 200*(1:k -1))-14
  #w_deseq2 <-deseq2_test(w, R$dat, 1:p, R$x, R$signal$signal)
  #save(w_deseq2, file=paste0(file.start, "_", n, "_deseq2.RData"))

  #n <- "w8b"
  #w <- cbind( 96 + 200*(1:k -1), 103 + 200*(1:k -1))
  #w_deseq2 <-deseq2_test(w, R$dat, 1:p, R$x, R$signal$signal)
  #save(w_deseq2, file=paste0(file.start, "_", n, "_deseq2.RData"))

}



run_win_tests2 <- function(file.start, p, waveQTL_loc, s0=c(0, 0, 0)){
  R <- getobj(paste0(file.start, "_fret.RData"))

  #Windows around peaks
  k <- p/200

  if(FALSE){
  n <- "w64o"
  w <- cbind( 69 + 200*(1:k -1), 132 + 200*(1:k -1))-14
  w_test <- window_test(w, dat=R$dat,
                            pos=1:p, x=R$x, signal=R$signal$signal, s0=s0)
  save(w_test, file=paste0(file.start, "_", n, ".RData"))
  w_waveqtl <-run_waveQTL(w, dat=R$dat, x=R$x, signal=R$signal$signal,
                              waveQTL_loc=waveQTL_loc)
  save(w_waveqtl, file=paste0(file.start, "_", n, "_wave.RData"))

  n <- "w8b"
  w <- cbind( 96 + 200*(1:k -1), 103 + 200*(1:k -1))
  w_test <- window_test(w, dat=R$dat,
                        pos=1:p, x=R$x, signal=R$signal$signal, s0=s0)
  save(w_test, file=paste0(file.start, "_", n, ".RData"))
  w_waveqtl <-run_waveQTL(w, dat=R$dat, x=R$x, signal=R$signal$signal,
                          waveQTL_loc=waveQTL_loc)
  save(w_waveqtl, file=paste0(file.start, "_", n, "_wave.RData"))
  }

  n <- "w32e"
  w <- cbind(seq(1, p - (p%%32), by=32), seq(32, p, by=32))
  w_test <- window_test(w, dat=R$dat,
                        pos=1:p, x=R$x, signal=R$signal$signal, s0=s0)
  save(w_test, file=paste0(file.start, "_", n, ".RData"))
  w_waveqtl <-run_waveQTL(w, dat=R$dat, x=R$x, signal=R$signal$signal,
                          waveQTL_loc=waveQTL_loc)
  save(w_waveqtl, file=paste0(file.start, "_", n, "_wave.RData"))


  n <- "w32b"
  w <- cbind( 85 + 200*(1:k -1), 116 + 200*(1:k -1))
  w_test <- window_test(w, dat=R$dat,
                        pos=1:p, x=R$x, signal=R$signal$signal, s0=s0)
  save(w_test, file=paste0(file.start, "_", n, ".RData"))
  w_waveqtl <-run_waveQTL(w, dat=R$dat, x=R$x, signal=R$signal$signal,
                          waveQTL_loc=waveQTL_loc)
  save(w_waveqtl, file=paste0(file.start, "_", n, "_wave.RData"))


}





