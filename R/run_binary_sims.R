#'@import dplyr

#'@export
run_bin <- function(seed, prefix, n, n.perm=500, waveQTL_loc="~/.local/bin/WaveQTL"){
  file.start <- paste0(prefix, "_", n)
  x <- rep(c(0, 1), each=15)
  g1 = function(x){return(list("ht"=5, "assoc"=0))}
  g2 <- function(x){return(list("ht"=rexp(n = 1, rate=1/5)+1.5, "assoc"=0))}
  g3 <- function(x){
    if(x==0) return( list("ht"=sample(c(1.5, 3, 7, 15), prob=c(0.55, 0.3, 0.1, 0.05), size=1), "assoc"=1))
    return(list("ht"=sample(c(1.5, 3, 7, 15), prob=c(0.2, 0.45, 0.25, 0.1), size=1), "assoc"=1))
  }
  pk.ht.funcs = c(g1, g2, g3)
  type.sequence=rep(c(1, 3, 1, 2, 3, 2), each=5)


  R <- cfdr_sims2(x, pk.ht.funcs, type.sequence, n.seg=c(2, 6),
                         seed=seed, n.perms=n.perm, s0=rep(0.05, 3),
                         level=c(0.02, 0.05, 0.1, 0.2),
                         save.data=TRUE, huber.maxit=50,
                         file.name=paste0(file.start, "_fret.RData"))
  R <- getobj(paste0(file.start, "_fret.RData"))

  #Even Windows
  w_50e <- cbind(seq(26, 6000-27, by=50), seq(75, 6000, by=50))
  w_50e_test <- window_test(w_50e, dat=R$dat, pos=1:6000,
                            x=R$x, signal=R$signal$signal, s0=c(0, 0, 0))
  save(w_50e_test, file=paste0(file.start, "_w50e.RData"))

  #Windows around peaks
  w_50b <- cbind( 76 + 200*(1:30 -1), 125 + 200*(1:30 -1))
  w_50b_test <- window_test(w_50b, dat=R$dat,pos=1:6000,
                            x=R$x, signal=R$signal$signal, s0=c(0, 0, 0))
  save(w_50b_test, file=paste0(file.start, "_w50b.RData"))

  #Width 64 bins for waveQTL
  #Even
  w_64e <- cbind(seq(1, 6000, by=64)[-94], seq(64, 6000, by=64))
  w_64e_test <- window_test(w_64e, dat=R$dat,
                            pos=1:6000, x=R$x, signal=R$signal$signal, s0=c(0, 0, 0))
  save(w_64e_test, file=paste0(file.start, "_w64e.RData"))
  w_64e_waveqtl <-run_waveQTL(w_64e, dat=R$dat, x=R$x, signal=R$signal$signal,
                      waveQTL_loc=waveQTL_loc)
  save(w_64e_waveqtl, file=paste0(file.start, "_w64e_wave.RData"))
  #At peaks
  w_64b <- cbind( 69 + 200*(1:30 -1), 132 + 200*(1:30 -1))
  w_64b_test <- window_test(w_64b, dat=R$dat,
                            pos=1:6000, x=R$x, signal=R$signal$signal, s0=c(0, 0, 0))
  save(w_64b_test, file=paste0(file.start, "_w64b.RData"))
  w_64b_waveqtl <-run_waveQTL(w_64b, dat=R$dat, x=R$x, signal=R$signal$signal,
                              waveQTL_loc=waveQTL_loc)
  save(w_64e_waveqtl, file=paste0(file.start, "_w64b_wave.RData"))

}

#'@export
run_bin2 <- function(seed, prefix, n, n.perm=500, waveQTL_loc="~/.local/bin/WaveQTL"){
  file.start <- paste0(prefix, "_", n)
  x <- rep(c(0, 1), each=15)
  g1 = function(x){return(list("ht"=5, "assoc"=0))}
  g2 <- function(x){return(list("ht"=rexp(n = 1, rate=1/5)+1.5, "assoc"=0))}
  g3 <- function(x){
    p1 = sample(c(1.5, 3, 7, 15), prob=c(0.55, 0.3, 0.1, 0.05), size=1)
    p2 = sample(c(1.5, 3, 7, 15), prob=c(0.2, 0.45, 0.25, 0.1), size=1)
    if(x==0) return( list("ht"=c(p1, p2), "assoc"=c(1, 1)))
    return(list("ht"=c(p2, p1), "assoc"=c(1, 1)))
  }
  pk.ht.funcs = c(g1, g2, g3)
  type.sequence=rep(c(1, 3, 1, 2, 3, 2), each=5)


  R <- cfdr_sims2(x, pk.ht.funcs, type.sequence, n.seg=c(2, 6),
                  seed=seed, n.perms=n.perm, s0=rep(0.05, 3),
                  level=c(0.02, 0.05, 0.1, 0.2),
                  save.data=TRUE, huber.maxit=50,
                  file.name=paste0(file.start, "_fret.RData"))
  R <- getobj(paste0(file.start, "_fret.RData"))

  #Even Windows
  w_50e <- cbind(seq(26, 6000-27, by=50), seq(75, 6000, by=50))
  w_50e_test <- window_test(w_50e, dat=R$dat, pos=1:6000,
                            x=R$x, signal=R$signal$signal, s0=c(0, 0, 0))
  save(w_50e_test, file=paste0(file.start, "_w50e.RData"))

  #Windows around peaks
  w_50b <- cbind( 76 + 200*(1:30 -1), 125 + 200*(1:30 -1))
  w_50b_test <- window_test(w_50b, dat=R$dat,pos=1:6000,
                            x=R$x, signal=R$signal$signal, s0=c(0, 0, 0))
  save(w_50b_test, file=paste0(file.start, "_w50b.RData"))

  #Width 64 bins for waveQTL
  #Even
  w_64e <- cbind(seq(1, 6000, by=64)[-94], seq(64, 6000, by=64))
  w_64e_test <- window_test(w_64e, dat=R$dat,
                            pos=1:6000, x=R$x, signal=R$signal$signal, s0=c(0, 0, 0))
  save(w_64e_test, file=paste0(file.start, "_w64e.RData"))
  w_64e_waveqtl <-run_waveQTL(w_64e, dat=R$dat, x=R$x, signal=R$signal$signal,
                              waveQTL_loc=waveQTL_loc)
  save(w_64e_waveqtl, file=paste0(file.start, "_w64e_wave.RData"))
  #At peaks
  w_64b <- cbind( 69 + 200*(1:30 -1), 132 + 200*(1:30 -1))
  w_64b_test <- window_test(w_64b, dat=R$dat,
                            pos=1:6000, x=R$x, signal=R$signal$signal, s0=c(0, 0, 0))
  save(w_64b_test, file=paste0(file.start, "_w64b.RData"))
  w_64b_waveqtl <-run_waveQTL(w_64b, dat=R$dat, x=R$x, signal=R$signal$signal,
                              waveQTL_loc=waveQTL_loc)
  save(w_64e_waveqtl, file=paste0(file.start, "_w64b_wave.RData"))
}
