
run_waveQTL <- function(windows,type.sequence, dat){

  n <- ncol(dat)
  p <- nrow(windows)
  S <- Intervals(windows)

  l <- rep(0, p)
  s <- get_signal(type.sequence, bandwidth=0)
  p0 <- nrow(s$signal)
  d <- distance_to_nearest(S, s$signal)
  l[d==0] <- 1
  N <- floor(runif(n=1, min=10, max=1e9))
  pvals <- c()
  for(i in 1:p){
    pheno.dat <- t(dat[windows[i, 1]:windows[i, 2], ])
    res <- WaveQTL_preprocess(Data = pheno.dat, meanR.thresh = 0)
    f <- tempfile(tmpdir = ".")
    write.table(res$WCs, file=paste0(f, "_pheno.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
    cat(res$filtered.WCs, file=paste0(f, "_use.txt"))
    cmd <- paste0("~/Desktop/Cluster_FDR/cfdrSims/WaveQTL-master/bin/WaveQTL -gmode 1 -g geno.txt -p ",
            f, "_pheno.txt -u ", f, "_use.txt -o temp", N, " -f ", n, " -numPerm 1000 -fph 2")
    system(cmd)
    pval <- read.table("output/temp.fph.pval.txt", header=TRUE)
    pvals <- c(pvals, pval[3, 1])
    unlink(paste0(f, "_pheno.txt"))
    unlink(paste0(f, "_use.txt"))
  }
  rates <- data.frame(t(sapply(sort(-log10(pvals)), FUN=function(x){
    tpr_nfp(s$signal, discoveries=windows[-log10(pvals) >= x, , drop=FALSE])
  })))
  return(list("pvals"=pvals, "rates"=rates))
}