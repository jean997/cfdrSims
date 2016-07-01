
#' Run Simulations
#'@description Run simulations
#'@param x Trait
#'@param pk.ht.funcs List of peak height functions. Each should take in a value of x
#'and return a a list of "height"=peak height (or vector of peak heights) and "assoc" as a
#'vector indicating if each peak is associated with x or not. "assoc" should be fixed for
#'all values of x.
#'@param type.sequence Vector of integers stating which peak height function to use.
#'@param seed Optional - set a seed.
#'@param n.perms Number of permutations to use to estimate lambda.
#'@param s0 Vector of length 3 giving additional variance to add to Poisson, Huber and t-test statistics.
#'@export
cfdr_sims2 <- function(x, pk.ht.funcs, type.sequence, n.seg=NULL,
                           seed=NULL, n.perms=500, s0=c(0, 0, 0),
                           level=c(0.02, 0.05, 0.1, 0.2),
                           save.data=FALSE, huber.maxit=50,
                            file.name=NULL){

  if(!is.null(seed)) set.seed(seed)
  stat.names = c("Poisson", "Huber", "t-test")

  if(all(x %in% c(0, 1))){
    pois_func = get_stats_pois_binary
    lm_func <- get_stats_ttest
  }else{
    pois_func = get_stats_pois_continuous
    lm_func=get_stats_lm
  }

  if(!is.null(file.name)){
    ff <- unlist(strsplit(file.name, ".RData"))[1]
    temp.name <- paste(ff, "_temp.RData", sep="")
  }

  stopifnot(all(sapply(pk.ht.funcs, FUN=class) == "function"))

  r <- length(type.sequence)
  p <- 200 * length(type.sequence)
  b <- length(level)

  #Segment Ends
  seg.ends=list()
  seg.ends[[1]] <- c(p)
  k <- 1
  if(!is.null(n.seg)){
    for(i in 1:length(n.seg)){
      stopifnot(p %% n.seg == 0)
      seg.ends[[i+1]] <- seq(p/n.seg[i], p, length.out=n.seg[i])
    }
    k <- k + length(n.seg)
  }

  #Build signal intervals object
  S <- cfdrSims:::get_signal2(pk.ht.funcs, type.sequence)

  perms <- replicate(n=n.perms, expr = {
    sample( x, size=length(x), replace=FALSE)
  })

  #n.f.disc <- n.t.disc <- tpr <-  array(0, dim=c(3, b, k))
  rates_at <- array(0, dim=c(3, k, b, 4))
  dimnames(rates_at)= list(stat.names, c(1, n.seg), level, c("tpr", "nfp", "ntp", "fdp"))
  #Generage Data
  P <- sapply(x, FUN=function(xx){
    ht.list = lapply(type.sequence, FUN=function(t){pk.ht.funcs[[t]](xx)$ht})
    yy <- unlist(lapply(ht.list, FUN=cfdrSims:::gen_profile))
    return(yy)
  })
  D <- apply(P, MARGIN=2, FUN=function(m){rpois(n=nrow(P), lambda=m)})

  if(!is.null(file.name)) save(D, file=temp.name)
  if(save.data){
    stats <- array(dim=c(3, p, 1+n.perms))
  }else{
    stats <- NULL
  }
  cl <- list()
  for(i in 1:length(stat.names)){
    cat(stat.names[i], "..")
    if(stat.names[i]=="Poisson"){
      Z <- pois_func(D, x, perms, s0=s0[i])
    }else if(stat.names[i]=="Huber"){
      Z <- get_stats_huber2(D, x, perms, s0=s0[i], maxit=huber.maxit)
    }else if(stat.names[i]=="t-test"){
      Z <- lm_func(D, x, perms, s0=s0[i])
    }

    if(save.data) stats[i, , ] <- Z
    if(!is.null(file.name)) save(D, stats, file=temp.name)
    cl[[i]] <- list()
    for(k in 1:length(seg.ends)){
      cl[[i]][[k]] <- get_clusters(Z, 1:p, bw=20, nlam=50,
                                   level=level, seg.ends=seg.ends[[k]])

      for(j in 1:b){
        if(nrow(cl[[i]][[k]]$clust[[j]])> 0){
          rates_at[i, k, j, ]<- cfdrSims:::tpr_nfp(Intervals(S$signal),
                                      discoveries=cl[[i]][[k]]$clust[[j]])
        }
      }
    }
    cat("\n")
  }

  if(!save.data) D <- NULL
  R <- list("rates_at"=rates_at,
            "stat.names"=stat.names, "cl"=cl, "signal"=S, "means"=P,
            "type.sequence"=type.sequence, "x"=x, "pk.pk.ht.funcs"=pk.ht.funcs,
            "dat"=D, "stats"=stats, "seed"=seed)
  if(!is.null(file.name)){
    save(R, file=file.name)
    unlink(temp.name)
  }else{
    return(R)
  }
}
