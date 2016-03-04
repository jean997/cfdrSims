
#' Run Simulations
#'@description Run simulations
#'@param type.sequence Vector of integers in 1:6 giving the type of each subregion.
#'@param B Number of simulations
#'@param sample.size Vector of length 2 giving sample sizes
#'@param seed Optional - set a seed.
#'@param n.perms Number of permutations to use to estimate lambda.
#'@param s0 Vector of length 3 giving additional variance to add to Poisson, Huber and t-test statistics.
#' @return A list with elements:
#' \describe{
#' \item{\code{n.false.discoveries}}{ Matrix 3 x B Giving total false discoveries for each method and each simulation.}
#' \item{\code{true.discoveries}}{ Array 3 x q x B where q is the number of regions of type 4:6}
#' }
#'@export
run_cfdr_sims <- function(type.sequence, B=10, sample.size=c(20, 20),
                     seed=NULL, n.perms=500, s0=c(0.1, 0.1, 0)){
  if(!is.null(seed)) set.seed(seed)

  q <- sum(type.sequence %in% 4:6)
  p <- 200 * length(type.sequence)
  cat(q, " regions with signal,", length(type.sequence)-q, " regions without.\n")
  if(q > 0){
    signal <- matrix(nrow=q, ncol=2)
    i <- 1; ct <- 1
    for(t in type.sequence){
      if(t %in% 4:6){
        signal[ct,] <- c(i+80, i + 118)
        ct <- ct + 1
      }
      i <- i + 200
    }
    signal <- Intervals(signal)
    not_signal <- interval_complement(signal)
  }else{
    not_signal <- Intervals(matrix( c(-Inf, Inf), ncol=2))
  }
  labs <- rep(c(0, 1), sample.size)
  perms <- replicate(n=n.perms, expr = {
    sample( labs, size=length(labs), replace=FALSE)
  })

  stat.names <- c("Poisson", "Huber", "t-test")
  n.false.discoveries <- matrix(nrow=3, ncol=B)
  if(q > 0) true.discoveries <- array(0, dim=c(3, q, B))
    else true.discoveries=NULL
  for(i in 1:B){
    cat(i, ": ")
    D <- sample_data(type.sequence=type.sequence, sample.size=sample.size)
    for(type in stat.names){
      cat(type, " ")
      j <- which(stat.names==type)
      if(type=="Poisson") Z <- get_stats_pois(D$dat , labs, perms, s0=s0[1])
        else if(type=="Huber") Z <- get_stats_huber(D$dat, labs, perms, s0=s0[2])
          else if(type=="t-test") Z <- get_stats_ttest(D$dat, labs, perms, s0=s0[3])
      cl <- get_clusters(Z, 1:p, 20, 20)
      if(cl$zsel < max(cl$z)){
        if(q > 0){
          td <- interval_overlap(cl$clust, signal)
          td_class <- unlist(lapply(td, FUN=length))
          td <- unlist(td)
          true.discoveries[j, td, i] <- 1
        }else{
          td_class <- rep(0, nrow(clpois$clust))
        }
        n.false.discoveries[j,i] <- sum(td_class==0)
      }else{
        n.false.discoveries[j, i] <- 0
      }
    }
    cat("\n")
  }
  return(list(n.false.discoveries, true.discoveries))
}
