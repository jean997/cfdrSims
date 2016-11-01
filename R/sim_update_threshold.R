
#'@export
sim_update_threshold <- function(prefix, n, ext="upd_rates", bandwidth=64,
                                 n.seg = c(), auto=c(50, 100, 150, 200, 300, 400)){
  file.name <- paste0(prefix, "_", n, "_fret.RData")
  save.name <- paste0(prefix, "_", n, "_", ext, ".RData")
  R <- getobj(file.name)
  p <- dim(R$dat)[1]
  n.perms <- dim(R$stats)[3]-1

  #segment.bounds
  K <- length(n.seg) + length(auto)
  sb=list()
  if(length(n.seg) > 0){
    for(i in 1:length(n.seg)){
      stopifnot(p %% n.seg == 0)
      sb[[i]] <- cbind(seq(1, p, by=p/n.seg[i]),  seq(p/n.seg[i], p, by=p/n.seg[i]))
    }
  }
  dd <- dim(R$rates)
  dd[2] <- K
  dn <- dimnames(R$rates)
  dn[[2]] <- c(as.character(n.seg), paste0("auto", auto))
  stat.names <- dimnames(R$rates)[[1]]
  level <- as.numeric(dimnames(R$rates)[[3]])
  b <- length(level)
  new.rates <- new.rates <- array(0, dim=dd)
  dimnames(new.rates) <- dimnames(new.rates) <- dn


  for(i in 1:length(stat.names)){
    cat(stat.names[i], " ")
    Z <- R$stats[i, , ]
    Zs <- apply(Z, MARGIN=2, FUN=function(y){
      ksmooth(x=1:p, y=y, x.points=1:p, bandwidth=64)$y
    })

    #Automatically determined intervals
    vv <- apply(Zs[,-1], MARGIN=1, FUN=var)
    for(k in 1:length(auto)){
      ss <- find_segments3(vv=vv, pos=1:p, min.length = auto[k],q=0.05)
      sb[[k + length(n.seg)]] <- data.frame("chrom"=rep("chr1", nrow(ss)),
                                            "start"=ss[,1], "stop"=ss[,2])
    }
    cat(" nseg: ")
    mtab <- R$mtabs[[i]]
    for(k in 1:K){
      cat(nrow(sb[[k]]), " ")
      fstep2 <- fret_step2(mtab$max1, mtab$max.perm, mtab$n.perm, mtab$zmin, sb[[k]])
      fstep3 <- fret_step3(fstep2, level)
      for(j in 1:b){
        if(level[j] %in% fstep3$Robs$fdr){
          ix <- which(fstep3$Robs$fdr == level[j])
          z <- as.numeric(rep(fstep3$z[1, ix, -c(1, 2)], fstep2$nbp))
          discoveries <- name_clusters_merged(x=Zs[,1], z=z, z0 = 0.3*mtab$zmin)
          new.rates[i, k, j, ]<- tpr_nfp(Intervals(R$signal$signal),
                                     discoveries=discoveries)
        }
      }
    }
    cat("\n")
  }
  #ret <- list("new.rates"=new.rates)
  save(new.rates, file=save.name)
  return(NULL)
}
