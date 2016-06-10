
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
#' \item{\code{n.f.disc}}{ Number of false disccoveries}
#' \item{\code{n.t.disc}}{ Array 3 x q x B where q is the number of regions of type 4:6}
#' }
#'@export
run_cfdr_sims1 <- function(type.sequence, sample.size=c(20, 20),
                     seed=NULL, n.perms=500, s0=c(0, 0, 0),
                     level=c(0.02, 0.05, 0.1, 0.2),
                     save.data=FALSE, type.def=NULL,
                     stat.names = c("Poisson", "Huber", "t-test")){

  if(!is.null(seed)) set.seed(seed)
  if(!is.null(type.def)) stopifnot(names(type.def) ==c("p1", "p2"))
    else type.def=define_types()

  r <- length(type.sequence)
  q <- sum(type.sequence %in% 4:6)
  p <- 200 * length(type.sequence)
  b <- length(level)
  cat(q, " regions with signal,", r-q, " regions without.\n")

  #Build signal intervals object
  S <- get_signal(type.sequence)

  labs <- rep(c(0, 1), sample.size)
  perms <- replicate(n=n.perms, expr = {
    sample( labs, size=length(labs), replace=FALSE)
  })

  n.f.disc <- n.t.disc <-  array(0, dim=c(3, b))

  D <- sample_data(type.sequence=type.sequence, sample.size=sample.size, type.def=type.def)
  if(save.data){
    stats <- array(dim=c(3, p, 1+n.perms))
  }else{
    stats <- NULL
  }
  for(i in 1:length(stat.names)){
    if(stat.names[i]=="Poisson") Z <- get_stats_pois(D$dat , labs, perms, s0=s0[i])
      else if(stat.names[i]=="Huber") Z <- get_stats_huber2(D$dat, labs, perms, s0=s0[i])
        else if(stat.names[i]=="t-test") Z <- get_stats_ttest(D$dat, labs, perms, s0=s0[i])

    cl <- get_clusters(Z, 1:p, bw=20, nlam=50, level=level)
    if(save.data) stats[i, , ] <- Z
    for(j in 1:b){
      if(nrow(cl$clust[[j]])> 0){
        rates <- cfdrSims:::tpr_nfp(S$signal, discoveries=cl$clust[[j]])
        n.t.disc[i, j] <- rates["ntp"]
        n.f.disc[i, j] <- rates["nfp"]
      }
    }
    cat("\n")
  }

  if(!save.data) D$dat <- NULL
  return(list("n.f.disc"=n.f.disc, "n.t.disc"=n.t.disc,
              "stat.names"=stat.names, "cl"=cl,
              "type.sequence"=type.sequence, "sample.size"=sample.size,
              "dat"=D$dat, "stats"=stats, "type.def"=type.def))
}
