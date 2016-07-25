
#All have t-test and box kernel smoother
#'@export
test_smooth_jade_sim <- function(data.file, profiles, bandwidth=20,
                                 smooth.type="box", thresh=1e-6){


  f <- getobj(data.file)
  labs <- rep(c(0, 1), f$sample.size)
  n <- nrow(f$Y)

  smooth.type <- match.arg(smooth.type, c("box", "spline"))
  if(smooth.type=="box"){
    smooth.func <- function(y){
      ksmooth(x=1:n, y=y, x.points = 1:n, bandwidth=bandwidth)$y
    }
  }else if(smooth.type=="spline"){
    smooth.func <- function(y){
      smooth.spline(x=1:n, y=y)$y
    }
  }

  #Build signal intervals object
  q0 <-rle( abs(profiles[,1]-profiles[,2]) > thresh )
  p0 <- length(q0$lengths)
  signal<- cbind(c(1, cumsum(q0$lengths)[-p0]+1)[q0$values], (cumsum(q0$lengths))[q0$values])
  signal <- Intervals(signal)
  #Point-wise signal
  pw_signal <- rep(0, nrow(profiles))
  pw_signal[ abs(profiles[,1]-profiles[,2]) > thresh] <- 1

  #Stats and smoothed stats
  stats <- cfdrSims:::t_stats(dat=f$Y, labs=labs, s0=0)
  ys <- smooth.func(stats)

  #Smoothed data and stats
  Ys <- apply(f$Y, MARGIN=2, FUN=function(y){
    smooth.func(y)
  })
  stat_Ys <- cfdrSims:::t_stats(dat=Ys, labs=labs, s0=0)

  #Rates for test then smooth
  yy <- sort(abs(ys), decreasing=TRUE)
  pw_ts <- data.frame(t(sapply(yy, FUN=function(x){
    sep <- as.numeric(abs(ys) > x)
    cbind(tpr.func(sep, pw_signal), fpr.func(sep, pw_signal))
  })))
  names(pw_ts) <- c("tpr", "fpr")
  yy <- yy[1:(0.5*n)]
  rw_ts <- data.frame(t(sapply(yy, FUN=function(x){
      cfdrSims:::tpr_nfp(signal, x=ys, z = x, z0=x)
      #z0=x means no merging
    })))
  rw_ts_merge <- data.frame(t(sapply(yy, FUN=function(x){
    cfdrSims:::tpr_nfp(signal, x=ys, z = x, z0=yy[length(yy)])
  })))


  #Rates for smooth then test
  yy <- sort(abs(stat_Ys), decreasing=TRUE)
  pw_st <- data.frame(t(sapply(yy, FUN=function(x){
    sep <- as.numeric(abs(stat_Ys) > x)
    cbind(tpr.func(sep, pw_signal), fpr.func(sep, pw_signal))
  })))
  names(pw_st) <- c("tpr", "fpr")
  yy <- yy[1:(0.5*n)]
  rw_st <- data.frame(t(sapply(yy, FUN=function(x){
    cfdrSims:::tpr_nfp(signal, x=stat_Ys, z = x, z0=x)
    #z0=x means no merging
  })))
  rw_st_merge <- data.frame(t(sapply(yy, FUN=function(x){
    cfdrSims:::tpr_nfp(signal, x=stat_Ys, z = x, z0=yy[length(yy)])
  })))


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
