
#All have t-test and box kernel smoother
#'@export
test_smooth_jade_sim <- function(file.prefix, which.rep,
                                 profiles, bandwidth=20,
                                 smooth.type=c("box", "spline"),
                                 stat.type=c("t", "huber"), s0=0,
                                 thresh=1e-6){

  stat.type <- match.arg(stat.type)
  #Build signal intervals object
  q0 <-rle( abs(profiles[,1]-profiles[,2]) > thresh )
  p0 <- length(q0$lengths)
  signal<- cbind(c(1, cumsum(q0$lengths)[-p0]+1)[q0$values], (cumsum(q0$lengths))[q0$values])
  signal <- Intervals(signal)
  #Point-wise signal
  pw_signal <- rep(0, nrow(profiles))
  pw_signal[ abs(profiles[,1]-profiles[,2]) > thresh] <- 1

  n <- nrow(profiles)


  smooth.type <- match.arg(smooth.type)
  if(smooth.type=="box"){
    smooth.func <- function(y){
      ksmooth(x=1:n, y=y, x.points = 1:n, bandwidth=bandwidth)$y
    }
  }else if(smooth.type=="spline"){
    smooth.func <- function(y){
      smooth.spline(x=1:n, y=y)$y
    }
  }

  if(stat.type=="t"){
    stat.func <- function(Y, labs){
      t_stats(Y, labs, s0=s0)
    }
  }else if(stat.type=="huber"){
    stat.func <- function(Y, labs){
      huber_stats2(Y, labs, s0=s0)
    }
  }
  #pw, rw, rw_merged
  ts_f <- list()
  ts_t <- list()
  ts_rates <- list()
  st_f <- list()
  st_t <- list()
  st_rates <- list()
  for(i in 1:3) ts_f[[i]] <- ts_t[[i]] <- st_f[[i]] <- st_t[[i]] <- list()
  j <- 1
  for(rep in which.rep){
    data.file <- paste0("data/", file.prefix, "_n", rep, "_data.RData")
    cat(data.file, "\n")
    f <- getobj(data.file)
    labs <- rep(c(0, 1), f$sample.size)
    #Stats and smoothed stats
    stats <-stat.func(f$Y, labs)
    ys <- smooth.func(stats)

    #Smoothed data and stats
    Ys <- apply(f$Y, MARGIN=2, FUN=function(y){
      smooth.func(y)
    })
    stat_Ys <- stat.func(Ys, labs)

    #Rates for test then smooth
    yy <- sort(abs(ys), decreasing=TRUE)
    pw_ts <- data.frame(t(sapply(yy, FUN=function(x){
      sep <- as.numeric(abs(ys) > x)
      cbind(tpr.func(sep, pw_signal), fpr.func(sep, pw_signal))
    })))
    names(pw_ts) <- c("tpr", "fpr")
    ts_t[[1]][[j]] <- pw_ts$tpr
    ts_f[[1]][[j]] <- pw_ts$fpr
    yy <- yy[1:(0.5*n)]
    rw_ts <- data.frame(t(sapply(yy, FUN=function(x){
      cfdrSims:::tpr_nfp(signal, x=ys, z = x)
      #z0=x means no merging
    })))
    ts_t[[2]][[j]] <- rw_ts$tpr
    ts_f[[2]][[j]] <- rw_ts$nfp
    rw_ts_merge <- data.frame(t(sapply(yy, FUN=function(x){
      cfdrSims:::tpr_nfp(signal, x=ys, z = x, z0=yy[length(yy)])
    })))
    ts_t[[3]][[j]] <- rw_ts_merge$tpr
    ts_f[[3]][[j]] <- rw_ts_merge$nfp

    #Rates for smooth then test
    yy <- sort(abs(stat_Ys), decreasing=TRUE)
    pw_st <- data.frame(t(sapply(yy, FUN=function(x){
      sep <- as.numeric(abs(stat_Ys) > x)
      cbind(tpr.func(sep, pw_signal), fpr.func(sep, pw_signal))
    })))
    names(pw_st) <- c("tpr", "fpr")
    st_t[[1]][[j]] <- pw_st$tpr
    st_f[[1]][[j]] <- pw_st$fpr
    yy <- yy[1:(0.5*n)]
    rw_st <- data.frame(t(sapply(yy, FUN=function(x){
      cfdrSims:::tpr_nfp(signal, x=stat_Ys, z = x, z0=x)
      #z0=x means no merging
    })))
    st_t[[2]][[j]] <- rw_st$tpr
    st_f[[2]][[j]] <- rw_st$nfp
    rw_st_merge <- data.frame(t(sapply(yy, FUN=function(x){
      cfdrSims:::tpr_nfp(signal, x=stat_Ys, z = x, z0=yy[length(yy)])
    })))
    st_t[[3]][[j]] <- rw_st_merge$tpr
    st_f[[3]][[j]] <- rw_st_merge$nfp
    j <- j + 1
  }

  st_rates[[1]] <- avg_by_interp(st_t[[1]], st_f[[1]])
  ts_rates[[1]] <- avg_by_interp(ts_t[[1]], ts_f[[1]])
  for(i in 2:3){
    st_rates[[i]] <- avg_by_ct(st_t[[i]], st_f[[i]])
    ts_rates[[i]] <- avg_by_ct(ts_t[[i]], ts_f[[i]])
  }


  R <- list("ts_t" = ts_t, "st_t"=st_t,
            "ts_f" = ts_f, "st_f"=st_f,
            "st_rates"=st_rates, "ts_rates" = ts_rates)
  return(R)
}

tpr.func <- function(x, labels){
  return(sum(x==1 & labels==1)/sum(labels==1))
}
fpr.func <- function(x, labels){
  return(sum(x==1 & labels==0)/sum(labels==0))
}



#'@export
test_smooth_fret_sim <- function(fret.file, bandwidth=20,
                                 sample.size=c(15, 15), stat.type=c("t", "huber"), s0=0,
                                 smooth.type=c("box", "spline")){

  stat.type <- match.arg(stat.type)
  smooth.type <- match.arg(smooth.type)
  R <- getobj(fret.file)
  #Build signal intervals object
  signal <- Intervals(R$signal$signal)
  #Point-wise signal
  pw_signal <- rep(0, nrow(R$dat))
  for(j in 1:nrow(signal)) pw_signal[signal[j, 1]:signal[j, 2]] <- 1

  n <- nrow(R$dat)

  if(smooth.type=="box"){
    smooth.func <- function(y){
      ksmooth(x=1:n, y=y, x.points = 1:n, bandwidth=bandwidth)$y
    }
  }else if(smooth.type=="spline"){
    smooth.func <- function(y){
      smooth.spline(x=1:n, y=y)$y
    }
  }

  if(stat.type=="t"){
    stat.func <- function(Y, labs){
      t_stats(Y, labs, s0=s0)
    }
  }else if(stat.type=="huber"){
    stat.func <- function(Y, labs){
      huber_stats2(Y, labs, s0=s0)
    }
  }

  labs <- rep(c(0, 1), sample.size)

  #Stats and smoothed stats
  stats <- stat.func(R$dat, labs)
  ys <- smooth.func(stats)

  #Smoothed data and stats
  Ys <- apply(R$dat, MARGIN=2, FUN=function(y){
    smooth.func(y)
  })
  stat_Ys <- stat.func(Ys, labs=labs)

  #Rates for test then smooth
  yy <- sort(abs(ys), decreasing=TRUE)
  pw_ts <- data.frame(t(sapply(yy, FUN=function(x){
    sep <- as.numeric(abs(ys) > x)
    cbind(tpr.func(sep, pw_signal), fpr.func(sep, pw_signal))
  })))
  names(pw_ts) <- c("tpr", "fpr")

  yy <- yy[1:(0.5*n)]
  rw_ts <- data.frame(t(sapply(yy, FUN=function(x){
    cfdrSims:::tpr_nfp(signal, x=ys, z = x)
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
  })))
  rw_st_merge <- data.frame(t(sapply(yy, FUN=function(x){
    cfdrSims:::tpr_nfp(signal, x=stat_Ys, z = x, z0=yy[length(yy)])
  })))

  results <- list("pw_ts"=pw_ts, "pw_st"  = pw_st,
                  "rw_ts"=rw_ts, "rw_st" = rw_st,
                  "rw_ts_merge" = rw_ts_merge, "rw_st_merge" = rw_st_merge)
  return(results)
}

#'@export
collect_test_smooth_fret <- function(file.start, which.rep, file.end){
  #pw, rw, rw_merged
  ts_f <- list()
  ts_t <- list()
  ts_rates <- list()
  st_f <- list()
  st_t <- list()
  st_rates <- list()
  for(i in 1:3) ts_f[[i]] <- ts_t[[i]] <- st_f[[i]] <- st_t[[i]] <- list()
  j <- 1
  for(rep in which.rep){
    R <- getobj(paste0(file.start, "_", rep, file.end))
    #Rates for test then smooth
    ts_t[[1]][[j]] <- R$pw_ts$tpr
    ts_f[[1]][[j]] <- R$pw_ts$fpr

    ts_t[[2]][[j]] <- R$rw_ts$tpr
    ts_f[[2]][[j]] <- R$rw_ts$nfp

    ts_t[[3]][[j]] <- R$rw_ts_merge$tpr
    ts_f[[3]][[j]] <- R$rw_ts_merge$nfp

    #Rates for smooth then test
    st_t[[1]][[j]] <- R$pw_st$tpr
    st_f[[1]][[j]] <- R$pw_st$fpr

    st_t[[2]][[j]] <- R$rw_st$tpr
    st_f[[2]][[j]] <- R$rw_st$nfp

    st_t[[3]][[j]] <- R$rw_st_merge$tpr
    st_f[[3]][[j]] <- R$rw_st_merge$nfp
    j <- j + 1
  }

  st_rates[[1]] <- avg_by_interp(st_t[[1]], st_f[[1]])
  ts_rates[[1]] <- avg_by_interp(ts_t[[1]], ts_f[[1]])
  for(i in 2:3){
    st_rates[[i]] <- avg_by_ct(st_t[[i]], st_f[[i]])
    ts_rates[[i]] <- avg_by_ct(ts_t[[i]], ts_f[[i]])
  }
  R <- list("ts_t" = ts_t, "st_t"=st_t,
            "ts_f" = ts_f, "st_f"=st_f,
            "st_rates"=st_rates, "ts_rates" = ts_rates)
  return(R)

}



