
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
cfdr_sims3 <- function(x, pk.ht.funcs, type.sequence,
                       seed=NULL, n.perms=500, s0=c(0, 0, 0),
                       n.seg=NULL, auto.min.length = NULL,
                       level=c(0.02, 0.05, 0.1, 0.2), huber.maxit=50,
                       save.data=FALSE, file.name=NULL,
                       random.peak.loc = FALSE, min.peak.sep=2){

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
  kk <- 1 + length(n.seg) + length(auto.min.length)

  #Segment bounds
  sb=list()
  sb[[1]] <- matrix(c(1, p), ncol=2)
  if(!is.null(n.seg)){
    for(i in 1:length(n.seg)){
      stopifnot(p %% n.seg == 0)
      sb[[i+1]] <- cbind(seq(1, p, by=p/n.seg[i]),  seq(p/n.seg[i], p, by=p/n.seg[i]))
    }
  }

  #Permutations
  perms <- replicate(n=n.perms, expr = {
    sample( x, size=length(x), replace=FALSE)
  })

  #Generage Data
  if(random.peak.loc){
    peak.starts <- sort(sample(1:p, size = length(type.sequence), replace=FALSE))
    sep <- 20 + min.peak.sep
    while(!all(diff(c(0, peak.starts, p)) >= sep )){
      ix <- min(min(which(diff(c(0, peak.starts, p)) < sep )), r)
      peak.starts[ix] <- sample((1:p)[-peak.starts], size=1)
      peak.starts <- sort(peak.starts)
    }
    type.sequence <- sample(type.sequence, size=r, replace=FALSE)
  }else{
    peak.starts <- seq(90, p-110, length.out=r)
  }
  P <- sapply(x, FUN=function(xx){
    peak.hts = sapply(type.sequence, FUN=function(t){pk.ht.funcs[[t]](xx)$ht})
    yy <- gen_profile2(peak.starts, peak.hts, total.length = p)
    return(yy)
  })
  D <- apply(P, MARGIN=2, FUN=function(m){rpois(n=nrow(P), lambda=m)})
  #Build signal object
  S <- get_signal3(pk.ht.funcs, type.sequence, peak.starts)

  #Output options
  rates <- array(0, dim=c(length(stat.names), kk, b, 4))
  dimnames(rates)= list(stat.names,
                        c(1, n.seg, paste0("auto", auto.min.length)),
                        level, c("tpr", "nfp", "ntp", "fdp"))
  if(!is.null(file.name)) save(D, file=temp.name)
  if(save.data){
    stats <- array(dim=c(3, p, 1+n.perms))
  }else{
    stats <- NULL
  }
  #cl will be  a list of length 3- one for each type of statistics
    #Each element of cl is a list of length kk (the number of segment partitions we consider)
    #each of this is an object produced by get_clusters2
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
    Zs <- apply(Z, MARGIN=2, FUN=function(y){
      ksmooth(x=1:p, y=y, x.points=1:p, bandwidth=20)$y
    })
    #Auto segment bounds
    vv <- apply(Zs[,-1], 1, var)
    if(!is.null(auto.min.length)){
      for(k in 1:length(auto.min.length)){
        sb[[k + length(n.seg) + 1]] <- find_segments(vv=vv, pos=1:p,
                                               min.length = auto.min.length[k],
                                               bandwidth=1, q=0.05)
      }
    }

    zmin <- as.numeric(quantile(Zs[,2:(n.perms + 1)], probs=c(0.95, 0.05)))
    z0 <- 0.3*zmin
    cl[[i]] <- list()
    for(k in 1:kk){
      cat(nrow(sb[[k]]), " ")
      cl[[i]][[k]] <- get_clusters2(Zs, pos=1:p, zmin=zmin, z0=z0,
                                level=c(0.02, 0.05, 0.1, 0.2), segment.bounds=sb[[k]])
      for(j in 1:b){
        if(nrow(cl[[i]][[k]]$clust[[j]])> 0){
          rates[i, k, j, ]<- tpr_nfp(Intervals(S$signal),
                                    discoveries=cl[[i]][[k]]$clust[[j]])
        }
      }
    }
    cat("\n")
  }

  if(!save.data) D <- NULL
  R <- list("rates"=rates,
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
