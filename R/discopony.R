


discopony_maxes1 <- function(dat.file, pheno.file, s0, z0, zmin,
                      strt, stp, seed, n.perm,
                      bandwidth=50, maxit=50,
                      chunk.size=10, prefix="",
                      bedops.loc="~/bedops/bin",
                      starch.loc="/projects/geneva/gcc-fs2/jean/Group_Curves/DNaseI/All_Download/jmorrison/SCLC-basal/perBase/"){

  stopifnot(length(chr.ends) ==2)


  X <- read_table(pheno.file, col_names=FALSE)
  n <- nrow(X)
  set.seed(seed)
  perms <- replicate(n=n.perm, expr = {
    sample( X[,2], size=n, replace=FALSE)
  })

#  strt0 <- strt <- chr.ends[1] + (which.chunk-1)*chunk.size*1000
#  stp <- strt <- chunk.size*1000
#  if(which.chunk > 1) strt0 <- strt - bandwidth
#  stp0 <- stp + bandwidth

#  file.name <- paste0(prefix, "chr", which.chr, "_chunk", which.chunk, ".RData")

#  system(paste0("./starch_extract.py --nostats --start ", strt0, " --stop ", stp0,
#        " --out ",  file.name, " --bedops-loc ", bdeops.loc,
#        "--starch-loc", starch.loc, " chr", which.chr[i],
#          " ", pheno.file))

  dat <- read_delim(dat.file, delim=" ")
  pos <- dat[,1]
  ix1 <- min(which(pos >=strt))
  ix2 <- max(which(pos <= stp))
  y <- huber_stats2(Y=dat[, -1], labs=X[,2],s0=s0, maxit=maxit)
  ys <- ksmooth_0(x=pos, y=y, bandwidth = bandwidth)[ix1:ix2]
  if(all(abs(ys) < zmin)){
    #unlink(file.name)
    return(0)
  }

  #Find peak heights
  #Intervals defined by z0:
  q0 <-rle( abs(ys) > z0 )
  p0 <- length(q0$lengths)
  ivls <- cbind(c(1, cumsum(q0$lengths)[-p0]+1)[q0$values], (cumsum(q0$lengths))[q0$values])

  max1 <- apply(ivls, MARGIN=1, FUN=function(iv){ max(abs(ys)[iv[1]:iv[2]])})

  perm.sm <- apply(perms, MARGIN=2, FUN=function(l){
      yy <- cfdrSims:::huber_stats2(dat[,-1], labs=l, s0=s0, maxit=maxit)
      cfdrSims:::ksmooth_0(x=pos, y=yy, bandwidth = bandwidth)[ix1:ix2]
  })
  max.perm <- apply(perm.sm, MARGIN=2, FUN=function(ys){
    if(all(abs(ys) <= z0)){
      c()
    }else{
      q0 <-rle( abs(ys) > z0 )
      p0 <- length(q0$lengths)
      ivls <- cbind(c(1, cumsum(q0$lengths)[-p0]+1)[q0$values], (cumsum(q0$lengths))[q0$values])
      apply(ivls, MARGIN=1, FUN=function(iv){ max(abs(ys)[iv[1]:iv[2]])})
    }
  })
  m <- sort(unlist(max.perm), decreasing=TRUE)
  mx <- cbind(m, (1:length(m))/(n.perm*(stp-strt + 1)))

  if(all(m < zmin)){
    mx <- cbind(zmin, 0)
  }else{
    mx <- mx[m >= zmin,]
  }
  file.name <- paste0(prefix, "chr", which.chr, "_chunk", which.chunk, "_mx.RData")
  id=paste0(which.chr,".", which.chunk)
  R <- list("max1"=max1, "mx"=mx, "id"=id)
  save(R, file=file.name)
}

discopony_choose_z <- function(file.list, zmin, nlam){
  log.lambda.min <- Inf
  log.lambda.max <- -Inf
  n.chunk <- 0
  for(f in file.list){
    R <- getobj(file.list)
    for(i in 1:length(R)){
      log.lambda.min <- min(log.lambda.min, log10(R[[i]]$mx[,2][mx[,2] > 0]))
      log.lambda.max <- max(log.lambda.max, log10(R[[i]]$mx[,2]))
    }
    n.chunk <- n.chunk + length(R)
  }
  lams <- seq(log.lambda.min, log.lambda.max, length.out=nlam)

  n.seg <- length(file.list)
  z <- matrix(nrow=nlam, ncol=n.seg+1)
  z[,1] <- lams
  for(f in file.list){

  }
  z <- data.frame(z)
  names(z) <- c("lambda", paste0("z", 1:nseg))

}
