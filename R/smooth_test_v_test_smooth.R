
#All have t-test and box kernel smoother
#'@export
test_smooth_jade_sim <- function(data.file, profiles, seed=NULL, bandwidth=20){


  if(!is.null(seed)) set.seed(seed)

  f <- getobj(data.file)
  labs <- rep(c(0, 1), f$sample.size)
  n <- nrow(f$Y)


  #Build signal intervals object
  q0 <-rle( abs(profiles[,1]-profiles[,2]) > 1e-6 )
  p0 <- length(q0$lengths)
  signal<- cbind(c(1, cumsum(q0$lengths)[-p0]+1)[q0$values], (cumsum(q0$lengths))[q0$values])
  signal <- Intervals(signal)
  #Point-wise signal
  pw_signal <- rep(0, nrow(profiles))
  pw_signal[ abs(profiles[,1]-profiles[,2]) > 1e-6 ] <- 1

  #Stats and smoothed stats
  stats <- t_stats(dat=f$Y, labs=labs, s0=0)
  ys <- ksmooth(x=1:n, y=stats, bandwidth=bandwidth, x.points=1:n)$y

  #Smoothed data and stats
  Ys <- apply(f$Y, MARGIN=2, FUN=function(y){
    ksmooth(x=1:n, y=y, bandwidth = bandwidth, x.points=1:n)$y
  })
  stat_Ys <- cfdrSims:::t_stats(dat=Ys, labs=labs, s0=0)

  #Rates for test then smooth
  yy <- sort(abs(ys), decreasing=TRUE)[1:(0.5*n)]
  rw_ts <- data.frame(t(sapply(yy, FUN=function(x){
      cfdrSims:::tpr_nfp(signal, x=ys, z = x, z0=x)
      #z0=x means no merging
    })))
  rw_ts_merge <- data.frame(t(sapply(yy, FUN=function(x){
    cfdrSims:::tpr_nfp(signal, x=ys, z = x, z0=yy[length(yy)])
  })))
  pw_ts <- data.frame(t(sapply(yy, FUN=function(x){
    sep <- as.numeric(abs(ys) > x)
    cbind(tpr.func(sep, pw_signal), fpr.func(sep, pw_signal))
  })))
  names(pw_ts) <- c("tpr", "fpr")

  #Rates for smooth then test
  yy <- sort(abs(stat_Ys), decreasing=TRUE)[1:(0.5*n)]
  rw_st <- data.frame(t(sapply(yy, FUN=function(x){
    cfdrSims:::tpr_nfp(signal, x=stat_Ys, z = x, z0=x)
    #z0=x means no merging
  })))
  rw_st_merge <- data.frame(t(sapply(yy, FUN=function(x){
    cfdrSims:::tpr_nfp(signal, x=stat_Ys, z = x, z0=yy[length(yy)])
  })))
  pw_st <- data.frame(t(sapply(yy, FUN=function(x){
    sep <- as.numeric(abs(stat_Ys) > x)
    cbind(tpr.func(sep, pw_signal), fpr.func(sep, pw_signal))
  })))
  names(pw_st) <- c("tpr", "fpr")



  R <- list("ys"=ys, "rw_ts"=rw_ts,
            "rw_ts_merge"=rw_ts_merge, "pw_ts"=pw_ts,
            "rw_st"=rw_st, "rw_st_merge"=rw_st_merge,
            "pw_st"=pw_st)
  return(R)
}

tpr.func <- function(x, labels){
  return(sum(x==1 & labels==1)/sum(labels==1))
}
fpr.func <- function(x, labels){
  return(sum(x==1 & labels==0)/sum(labels==0))
}
