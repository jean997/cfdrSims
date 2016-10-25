#'@import wavethresh


#'@export
run_waveQTL <- function(windows, dat, x, signal, level=c(0.02, 0.05, 0.1, 0.2),
                        waveQTL_loc="~/Desktop/Cluster_FDR/cfdrSims/WaveQTL-master/bin/WaveQTL"){

  N <- floor(runif(n=1, min=10, max=1e9)) #Random number label
  geno <- c("chr1.1", "A", "G", x) #''genotype'' info
  cat(geno, file=paste0("geno_", N, ".txt"))

  n <- ncol(dat)
  p <- nrow(windows)
  w<- Intervals(windows)
  if(!class(signal)=="Intervals") signal <- Intervals(signal)

  #Window labels (signal or not)
  l <- rep(0, p)
  d <- distance_to_nearest(w, signal)
  l[d==0] <- 1


  pvals <- c()
  for(i in 1:p){
    pheno.dat <- t(dat[windows[i, 1]:windows[i, 2], ])
    res <- WaveQTL_preprocess(Data = pheno.dat, meanR.thresh = 0)
    f <- tempfile(tmpdir = ".")
    write.table(res$WCs, file=paste0(f, "_pheno.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
    cat(res$filtered.WCs, file=paste0(f, "_use.txt"))
    cmd <- paste0(waveQTL_loc, " -gmode 1 -g geno_", N, ".txt -p ",
            f, "_pheno.txt -u ", f, "_use.txt -o temp", N, " -f ", n, " -numPerm 100000 -fph 2")
    system(cmd)
    pval <- read.table(paste0("output/temp", N, ".fph.pval.txt"), header=TRUE)
    pvals <- c(pvals, pval[3, 1])
    unlink(paste0(f, "_pheno.txt"))
    unlink(paste0(f, "_use.txt"))
  }
  qvals <- p.adjust(pvals, method="BH")
  rates <- data.frame(t(sapply(sort(-log10(pvals)), FUN=function(xx){
    tpr_nfp(signal, discoveries=w[-log10(pvals) >= xx, , drop=FALSE])
  })))
  rates_at <- t(sapply(level, FUN=function(ll){
    tpr_nfp(signal=signal, discoveries = w[qvals<= ll, ])
  }))

  unlink(paste0("geno_", N, ".txt"))
  return(list("pvals"=pvals, "qvals"=qvals, "rates"=rates, "rates_at"=rates_at))
}
