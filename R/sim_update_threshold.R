
#'@export
sim_update_threshold <- function(prefix, n, ext="upd_rates.RData",
                                 n.seg = c(), auto=c(50, 100, 150, 200, 300, 400),
                                 vv.bandwidth=1){
  file.name <- paste0(prefix, "_", n, "_fret.RData")
  save.name <- paste0(prefix, "_", n, "_", ext)
  R <- getobj(file.name)
  p <- dim(R$stats)[2]
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
  new.rates.sgn <- new.rates.usgn <- array(0, dim=dd)
  dimnames(new.rates.sgn) <- dimnames(new.rates.usgn) <- dn

  cl_unsigned <- list()
  cl_signed <- list()
  for(i in 1:length(stat.names)){
    cat(stat.names[i], " ")
    Z <- R$stats[i, , ]
    Zs <- apply(Z, MARGIN=2, FUN=function(y){
      ksmooth(x=1:p, y=y, x.points=1:p, bandwidth=20)$y
    })

    zmin_sgn <- as.numeric(quantile(Zs[,-1], probs=c(0.95, 0.05)))
    z0_sgn <- 0.3*zmin_sgn

    zmin_usgn <- as.numeric(quantile(abs(Zs[,-1]), probs=0.9))
    z0_usgn <- 0.3*zmin_usgn

    #Automatically determined intervals
    vv <- apply(Zs[,-1], MARGIN=1, FUN=var)
    for(k in 1:length(auto)){
      sb[[k + length(n.seg)]] <- find_segments(vv=vv, pos=1:p, min.length = auto[k],
                                               bandwidth=vv.bandwidth, q=0.05)
    }
    cl_signed[[i]] <- cl_unsigned[[i]] <- list()
    cat(" nseg: ")
    for(k in 1:K){
      cat(nrow(sb[[k]]), " ")
      cl_signed[[i]][[k]] <- get_clusters2(Zs, 1:p, zmin_sgn, z0_sgn,
                                    level=level, segment.bounds=sb[[k]])
      cl_unsigned[[i]][[k]] <- get_clusters2(Zs, 1:p, zmin_usgn, z0_usgn,
                                           level=level, segment.bounds=sb[[k]])
      for(j in 1:b){
        if(nrow(cl_signed[[i]][[k]]$clust[[j]])> 0){
          new.rates.sgn[i, k, j, ]<- tpr_nfp(Intervals(R$signal$signal),
                                                discoveries=cl_signed[[i]][[k]]$clust[[j]])
        }
        if(nrow(cl_unsigned[[i]][[k]]$clust[[j]])> 0){
          new.rates.usgn[i, k, j, ]<- tpr_nfp(Intervals(R$signal$signal),
                                     discoveries=cl_unsigned[[i]][[k]]$clust[[j]])
        }
      }
    }
    cat("\n")
  }
  ret <- list("cl_signed"=cl_signed, "cl_unsigned"=cl_unsigned,
              "new.rates.usgn"=new.rates.usgn, "new.rates.sgn"=new.rates.sgn )
  save(ret, file=save.name)
  return(NULL)
}
