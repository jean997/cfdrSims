

rates_from_jade_sim <- function(data.file, profiles, n.perm=0, seed=NULL, huber.maxit=50,
                                bandwidth=20, level=0.1,
                                stat.names = c("Poisson", "Huber", "t-test"),
                                s0=rep(0, length(stat.names))){

  if(!is.null(seed)) set.seed(seed)


  f <- getobj(data.file)
  labs <- rep(c(0, 1), f$sample.size)
  n <- nrow(f$Y)

  if(n.perm > 0){
    perms <- replicate(n=n.perm, expr = {
      sample( labs, size=length(labs), replace=FALSE)
    })
  }else{
    perms=NULL
  }


  #Build signal intervals object
  q0 <-rle( abs(profiles[,1]-profiles[,2]) > 1e-6 )
  p0 <- length(q0$lengths)
  signal<- cbind(c(1, cumsum(q0$lengths)[-p0]+1)[q0$values], (cumsum(q0$lengths))[q0$values])
  signal <- Intervals(signal)

  stats <- array(dim=c(3, n, 1+n.perm))
  rates <- array(dim=c(3, n, 4))
  cl <- list()
  for(i in 1:length(stat.names)){
    if(stat.names[i]=="Poisson") stats[i, , ] <- get_stats_pois(f$Y , labs, perms, s0=s0[i])
    else if(stat.names[i]=="Huber") stats[i, , ] <- get_stats_huber2(f$Y, labs, perms, s0=s0[i], maxit=huber.maxit)
    else if(stat.names[i]=="t-test") stats[i, , ] <- get_stats_ttest(f$Y, labs, perms, s0=s0[i])

    cl[[i]] <- get_clusters(stats[i, , ], 1:n, bw=bandwidth, nlam=50,
                            level=level, z.min.quantile = 0.8)

    ys <- ksmooth(x=1:n, y=stats[i, , 1], bandwidth=bandwidth, x.points=1:n)$y

    rates[i, , ] <- t(sapply(cl[[i]]$z[, 2], FUN=function(x){
      cfdrSims:::tpr_nfp(signal, x=ys, z = x, z0=cl[[i]]$z0)
    }))
  }

  R <- list("ys"=ys, "rates"=rates,
            "stat.names"=stat.names, "cl"=cl,
            "file"=dat.file, "stats"=stats)
  return(R)

}
