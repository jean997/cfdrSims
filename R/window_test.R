
#' Test simulations in windows
#'@description Calculate test statistics binning over aupplied regions
#'@param windows q by 2 matrix giving region boundaries
#'@param type.sequence Vector of integers in 1:6 giving the type of each subregion.
#'@param dat p x n matrix of data
#'@param labs 0/1 treatment status of each sample (length n)
#'@param s0
#'@param stat.names Length k. Contains elements of "Poisson", "Huber", "t-test"
#' @return A list with elements:
#' \describe{
#' \item{\code{stats}}{ matrix k by q giving statistic for each bin }
#' \item{\code{rate_list}}{list of length k. Each element is a three column matrix with columns fdp, tp, fp}
#' }
#'@export
window_test <- function(windows, type.sequence, dat, labs, s0=c(0, 0, 0),
                                stat.names = c("Poisson", "Huber", "t-test")){

  n <- ncol(dat)
  p <- nrow(windows)
  S <- Intervals(windows)

  l <- rep(0, p)
  s <- cfdrSims:::get_signal(type.sequence, bandwidth=0)
  p0 <- nrow(s$signal)
  d <- distance_to_nearest(S, s$signal)
  l[d==0] <- 1

  window.dat <- matrix(nrow=p, ncol=n)
  stats <- qvals <- pvals <- matrix(nrow=length(stat.names), ncol=p)
  for(i in 1:p){
    window.dat[i, ] <- colSums(dat[windows[i, 1]:windows[i, 2], ])
  }
  rate_list <- list()
  for(i in 1:length(stat.names)){
    if(stat.names[i]=="Poisson") stats[i,] <- pois_regression(Y=window.dat, labs=labs, s0=s0[i])
    else if(stat.names[i]=="Huber")  stats[i, ] <- cfdrSims:::huber_stats2(Y=window.dat, labs=labs, s0=s0[i])
    else if(stat.names[i]=="t-test") stats[i, ]<- cfdrSims:::t_stats(window.dat, labs=labs, s0=s0[i])
    rate_list[[i]] <- data.frame(t(sapply(sort(abs(stats[i,])), FUN=function(x){
      tpr_nfp(s$signal, discoveries=S[abs(stats[i,]) >= x, , drop=FALSE])
    })))
    pvals[i,] <- sapply(stats[i,], FUN=cfdrSims:::get.pt, df=n-2)
    qvals[i,] <- p.adjust(pvals[i,], method="BH")
  }
  return(list("stats"=stats, "rate_list"=rate_list, "stat.names"=stat.names, "pvals"=pvals, "qvals"=qvals))



}
