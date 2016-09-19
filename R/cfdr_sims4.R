
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
cfdr_sims4 <- function(x, pk.ht.funcs, type.sequence,
                       seed=NULL, n.perms=500,
                       n.seg=c(), auto.min.length = c(),
                       level=c(0.02, 0.05, 0.1, 0.2),
                       save.data=FALSE, file.name=NULL, data.only=FALSE,
                       random.peak.loc = FALSE, min.peak.sep=2, s0=c(0, 0, 0),
                       stat.funcs = c(qpois_stats_binary, huber_stats, t_stats),
                       stat.names=c("Poisson", "Huber", "T")){

  if(!is.null(seed)) set.seed(seed)

  if(!all(x %in% c(0, 1))){
    if(qpois_stats_binary %in% stats.funcs) stop("qpois_stats_binary onlly for binary traits.\n")
    if(t_stats %in% stats.funcs) stop("t_stats onlly for binary traits.\n")
  }
  stopifnot(length(s0)==length(stat.funcs))

  if(!is.null(file.name)){
    ff <- unlist(strsplit(file.name, ".RData"))[1]
    temp.name <- paste(ff, "_temp.RData", sep="")
  }

  stopifnot(all(sapply(pk.ht.funcs, FUN=class) == "function"))
  stopifnot(all(sapply(stat.funcs, FUN=class) == "function"))
  stopifnot(length(stat.funcs)==length(stat.names))
  r <- length(type.sequence)
  p <- 200 * length(type.sequence)
  b <- length(level)
  kk <- 1 + length(n.seg) + length(auto.min.length)

  #Segment bounds
  sb=list()
  sbnames <- c("1")
  sb[[1]] <- data.frame("chrom"="chr1", "start"=1, "stop"=p)
  if(length(n.seg)> 0){
    for(i in 1:length(n.seg)){
      stopifnot(p %% n.seg == 0)
      sb[[i+1]] <- data.frame("chrom"=rep("chr1", n.seg[i]),
                              "start"=seq(1, p, by=p/n.seg[i]),
                              "stop"=seq(p/n.seg[i], p, by=p/n.seg[i]))
    }
    sbnames <- c(sbnames, as.character(n.seg))
  }
  if(length(auto.min.length) > 0) sbnames <- c(sbnames, paste0("auto", auto.min.length))

  #Permutations
  perms <- replicate(n=n.perms, expr = {
    sample( x, size=length(x), replace=FALSE)
  })
  perms <- cbind(x, perms)

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
    yy <- cfdrSims:::gen_profile2(peak.starts, peak.hts, total.length = p)
    return(yy)
  })
  D <- apply(P, MARGIN=2, FUN=function(m){rpois(n=nrow(P), lambda=m)})
  #Build signal object
  S <- get_signal3(pk.ht.funcs, type.sequence, peak.starts)
  if(data.only){
    ret=list("dat"=D, "means"=P, "signa"=S, "x"=x)
    save(ret, file=file.name)
    return(ret)
  }
  #Output options
  rates <- array(0, dim=c(length(stat.funcs), kk, b, 4))
  dimnames(rates)= list(stat.names, sbnames,
                        level, c("tpr", "nfp", "ntp", "fdp"))
  if(!is.null(file.name)) save(D, file=temp.name)

  if(save.data){
    stats <- array(dim=c(3, p, n.perms+1))
  }else{
    stats <- NULL
  }


  #mtabs will be  a list of length 3- one for each type of statistics
  mtabs <- list()
  for(i in 1:length(stat.names)){
    cat(stat.names[i], "..")
    Z <- apply(perms, MARGIN=2, FUN=function(l){
      stat.funcs[[i]](D, l, s0=s0[i])[3,]
    })
    if(save.data) stats[i, , ] <- Z
    if(!is.null(file.name)) save(D, stats, file=temp.name)
    Zs <- apply(Z, MARGIN=2, FUN=function(y){
      ksmooth(x=1:p, y=y, x.points=1:p, bandwidth=20)$y
    })

    #zmin <- as.numeric(quantile(Zs[,-1], probs=c(0.95, 0.05)))
    zmin <- as.numeric(quantile(abs(Zs[,-1]), probs=c(0.9)))
    #Maxes table

    mtab <- maxes_table_for_sims(Zs, 1:p, zmin)
    mtabs[[i]] <- mtab
    #Auto segment bounds
    vv <- apply(Zs[,-1], 1, var)
    if(length(auto.min.length)> 0){
      for(k in 1:length(auto.min.length)){
        ss <- find_segments(vv=vv, pos=1:p, min.length = auto.min.length[k],
                            bandwidth=1, q=0.05)
        sb[[k + length(n.seg) + 1]]  <- data.frame("chrom"=rep("chr1", nrow(ss)),
                                                   "start"=ss[,1], "stop"=ss[,2])
      }
    }

    for(k in 1:kk){
      cat(nrow(sb[[k]]), " ")
      fstep2 <- fret_step2(mtab$max1, mtab$max.perm, mtab$n.perm, zmin, sb[[k]])
      fstep3 <- fret_step3(fstep2, level)
      for(j in 1:b){
        if(level[j] %in% fstep3$Robs$fdr){
          ix <- which(fstep3$Robs$fdr == level[j])
          z <- as.numeric(rep(fstep3$z[1, ix, -c(1, 2)], fstep2$nbp))
          discoveries <- name_clusters_merged(x=Zs[,1], z=z, z0 = 0.3*zmin)
          rates[i, k, j, ]<- tpr_nfp(Intervals(S$signal),
                                     discoveries=discoveries)
        }
      }
    }
    cat("\n")
  }

  if(!save.data) D <- NULL
  R <- list("rates"=rates,
            "stat.names"=stat.names, "mtabs"=mtabs, "signal"=S, "means"=P,
            "type.sequence"=type.sequence, "x"=x, "pk.pk.ht.funcs"=pk.ht.funcs,
            "dat"=D, "stats"=stats, "seed"=seed)
  if(!is.null(file.name)){
    save(R, file=file.name)
    unlink(temp.name)
  }else{
    return(R)
  }
}
