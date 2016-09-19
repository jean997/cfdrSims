
#' Test simulations in windows
#'@description Calculate test statistics binning over aupplied regions
#'@param windows q by 2 matrix giving region boundaries
#'@param pos Positions
#'@param dat p x n matrix of data
#'@param x treatment status of each sample (length n)
#'@param signal
#'@param s0
#'@param level
#' @return A list
#'@export
window_test <- function(windows, dat, pos, x, signal,
                        level=c(0.02, 0.05, 0.1, 0.2),
                        s0=c(0, 0, 0),
                        stat.funcs = c(qpois_stats_binary, huber_stats, t_stats),
                        stat.names=c("Poisson", "Huber", "T")){

  stopifnot(length(s0)==length(stat.funcs))
  stopifnot(all(sapply(stat.funcs, FUN=class) == "function"))
  stopifnot(length(stat.funcs)==length(stat.names))

  n <- ncol(dat)
  p <- nrow(windows)
  w<- Intervals(windows)
  if(!class(signal)=="Intervals") signal <- Intervals(signal)

  #Window labels (signal or not)
  l <- rep(0, p)
  d <- distance_to_nearest(w, signal)
  l[d==0] <- 1

  #Window data
  win <- sapply(pos, FUN=function(x){
      w1 <- which(x >= windows[,1])
      w2 <- which(x <= windows[,2])
      wov <- intersect(w1, w2)
      if(length(wov)==0) return(0)
        else if (length(wov)==1) return(wov)
          else return(NA)
  })
  dat <- data.frame(dat)
  dat$win <- win
  window.dat <- dat %>% group_by(win) %>% summarise_each(funs(sum))
  window.dat <- window.dat[window.dat$win > 0,]

  stats <- qvals <- pvals <- matrix(nrow=length(stat.names), ncol=p)
  rates <- array(dim=c(length(stat.names), p, 4))
  rates_at <- array(dim=c(length(stat.names), length(level), 4 ))
  dimnames(rates) = list(stat.names, 1:p, c("tpr", "nfp", "ntp", "fdp"))
  dimnames(rates_at)= list(stat.names, level, c("tpr", "nfp", "ntp", "fdp"))
  for(i in 1:length(stat.names)){

    stats[i,] <- stat.funcs[[i]](window.dat[,-1], x, s0=s0[i])[3,]

    rates[i, , ] <- t(sapply(sort(abs(stats[i,])), FUN=function(xx){
      tpr_nfp(signal, discoveries=w[abs(stats[i,]) >= xx, , drop=FALSE])
    }))
    pvals[i,] <- sapply(stats[i,], FUN=function(ss){
      2*pt(abs(ss), df=n-2, lower.tail=FALSE)
    })
    qvals[i,] <- p.adjust(pvals[i,], method="BH")

    rates_at[i, , ] <- t(sapply(level, FUN=function(ll){
      tpr_nfp(signal=signal, discoveries = w[qvals[i, ] <= ll, ])
    }))
  }
  return(list("stats"=stats, "rates"=rates, "rates_at"=rates_at, "windows"=windows,
              "stat.names"=stat.names, "pvals"=pvals, "qvals"=qvals, "class"=l))

}

#'@export
deseq2_test <- function(windows, dat, pos, x, signal,
                      level=c(0.02, 0.05, 0.1, 0.2)){


  n <- ncol(dat)
  p <- nrow(windows)
  w<- Intervals(windows)
  if(!class(signal)=="Intervals") signal <- Intervals(signal)

  #Window labels (signal or not)
  l <- rep(0, p)
  d <- distance_to_nearest(w, signal)
  l[d==0] <- 1

  #Window data
  win <- sapply(pos, FUN=function(x){
    w1 <- which(x >= windows[,1])
    w2 <- which(x <= windows[,2])
    wov <- intersect(w1, w2)
    if(length(wov)==0) return(0)
    else if (length(wov)==1) return(wov)
    else return(NA)
  })
  dat <- data.frame(dat)
  dat$win <- win
  window.dat <- dat %>% group_by(win) %>% summarise_each(funs(sum))
  window.dat <- window.dat[window.dat$win > 0,]

  counts <- as.matrix(window.dat[,-1])
  pheno <- data.frame(cbind(1:length(x), x))
  if(all(x%in% c(0, 1))) pheno$x <- factor(pheno$x)
  names(pheno) <- c("name", "x")
  colnames(counts) = pheno$name
  deseq_obj <- DESeqDataSetFromMatrix(counts, pheno, as.formula("~x"))
  norm_table=matrix(1, nrow=nrow(counts), ncol=ncol(counts))
  normalizationFactors(deseq_obj) = norm_table
  deseq_obj = DESeq(deseq_obj)
  results = results(deseq_obj, independentFiltering =FALSE)

  rates_at <- t(sapply(level, FUN=function(ll){
    tpr_nfp(signal=signal, discoveries = w[results$padj <= ll, ])
  }))
  rates <- t(sapply(sort(abs(results$stat)), FUN=function(xx){
    tpr_nfp(signal, discoveries=w[abs(results$stat) >= xx, , drop=FALSE])
  }))
  return(list("results"=results, "rates"=rates, "rates_at"=rates_at, "pvalue" = results$pvalue,
              "windows"=windows, "class"=l))
}
