


discopony_maxes <- function(pheno.file, chr.ends, s0, z0,
                      seed, n.perm, bandwidth=20, maxit=50,
                      which.chr=1:nrow(chr.ends), nlam=50,
                      chunk.size=2){

  stopifnot(length(which.chr)==nrow(chr.ends))

  X <- read_table(pheno.file, col_names=FALSE)

  set.seed(seed)
  perms <- replicate(n=n.perm, expr = {
    sample( X[,2], size=length(labs), replace=FALSE)
  })

  n.chrom <- nrow(chr.ends)

  myfiles <- c()
  max.ll <- 0
  min.ll <- Inf
  for(i in 1:n.chrom){
    file.start <- tempfile(tmpdir=".", pattern=paste0("chr", which.chr[i]))
    n.chunk <- ceiling((chr.ends[i,2]-chr.ends[i, 1])/(chunk.size*1000))
    strt1 <- seq(chr.ends[i, 1], chr.ends[i, 2], by=chunk.size*1000)
    for(j in 1:n.chunk){
      file.name <- paste0(file.start, "_chunk", j, ".RData")
      if(j == 1){
        strt <- strt1[j]
        stp <- strt1[j] + chunk.size*1000 + bandwidth
        x.range=c(1, chunk.size*1000)
      }else{
        strt <- strt1[j] - bandwidth
        stp <- strt1[j] + (chunk.size*1000) + bandwidth
        x.range <- c(bandwidth+1, bandwidth+1+(chunk.size*1000) )
      }
      myfiles[j] = tempfile()
      system(paste0("./starch_extract.py --start ", strt, " --stop ", stp,
        " --out ",  file.name, " --bedops-loc ~/bedops/bin/ chr", which.chr[i],
          " ", pheno.file))

      dat <- read_table(myfiles[j])
      Z <- get_stats_huber2(dat[,-1], X[,2], perms, s0=s0, maxit=maxit)
      xs <- ksmooth(x=dat[,1], y=Z[,1], bandwidth = bandwidth,
                    x.points=dat[x.range[1]:x.range[2], 1])$y
      if(all(abs(xs) < z0)){
        unlink(file.name)
        next
      }
      mx <- choose_z_even(perm.stats=Z[,-1], nlam=nlam,
                                           bw=bandwidth, pos=dat[,1], z0=z0,
                                          x.range=x.range, return.z=FALSE,
                                          keep.lists=TRUE)[[1]]
      #min.ll <- max(min.ll, min(log10(mx[,2])))
      #max.ll <- min(max.ll, max(log10(mx[,2])))
      save(mx, file=file.name)
    }
  }

}
