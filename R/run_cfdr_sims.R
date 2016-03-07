
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
                     seed=NULL, n.perms=500, s0=c(0.1, 0.1, 0),
                     save.data=FALSE, type.def=NULL,
                     stat.names = c("Poisson", "Huber", "t-test")){
  if(!is.null(seed)) set.seed(seed)
  if(!is.null(type.def)) stopifnot(names(type.def) ==c("p1", "p2"))

  r <- length(type.sequence)
  q <- sum(type.sequence %in% 4:6)
  p <- 200 * length(type.sequence)
  cat(q, " regions with signal,", r-q, " regions without.\n")
  #Build signal intervals object
  if(q > 0){
    signal <- matrix(nrow=q, ncol=2)
    i <- 1; ct <- 1
    for(t in type.sequence){
      if(t %in% 4:6){
        signal[ct,] <- c(i+80, i + 118) #Signal is spread out due to smoothing
        ct <- ct + 1
      }
      i <- i + 200
    }
    signal <- Intervals(signal)
    not_signal <- interval_complement(signal)
  }else{
    not_signal <- Intervals(matrix( c(-Inf, Inf), ncol=2))
  }
  #Intervals object giving each region
  rI <- Intervals(cbind(200*((1:r)-1) + 1, 200*(1:r)))

  labs <- rep(c(0, 1), sample.size)
  perms <- replicate(n=n.perms, expr = {
    sample( labs, size=length(labs), replace=FALSE)
  })

  n.f.disc <- array(dim=c(3, r, B))
  if(q > 0) n.t.disc <- array(0, dim=c(3, q, B))
    else n.t.disc=NULL
  if(save.data){
    dat <- array(dim=c(p, sum(sample.size), B))
    stats <- array(dim=c(3, p, 1+n.perms, B))
  }else{
    dat <- NULL
    stats <- NULL
  }
  ref.names <- c("Poisson", "Huber", "t-test")
  for(i in 1:B){
    cat(i, ": ")
    D <- sample_data(type.sequence=type.sequence, sample.size=sample.size, type.def=type.def)
    if(save.data) dat[,,i] <- D$dat
    for(st in stat.names){
      cat(st, " ")
      j <- which(ref.names==st)
      if(st=="Poisson") Z <- get_stats_pois(D$dat , labs, perms, s0=s0[1])
        else if(st=="Huber") Z <- get_stats_huber(D$dat, labs, perms, s0=s0[2])
          else if(st=="t-test") Z <- get_stats_ttest(D$dat, labs, perms, s0=s0[3])
      cl <- get_clusters(Z, 1:p, bw=20, nz=20)
      if(save.data) stats[j, , , i] <- Z
      if(cl$zsel < max(cl$z)){
        if(q > 0){
          td <- interval_overlap(cl$clust, signal)
          #1 if overlapping with signal, 0 otherwise
          td_class <- unlist(lapply(td, FUN=length))
          td <- unlist(td)
          n.t.disc[j, td, i] <- 1
        }else{
          td_class <- rep(0, nrow(cl$clust))
        }
        td_region <- unlist(interval_overlap(cl$clust, rI))
        f.disc.region <- rep(0, r)
        f.disc.region[td_region[td_class==0]] <- 1
        n.f.disc[j, ,i] <- f.disc.region
      }else{
        n.f.disc[j, ,i] <- rep(0, r)
      }
    }
    cat("\n")
  }
  if(is.null(type.def)) type.def=define_types()
  return(list("n.f.disc"=n.f.disc, "n.t.disc"=n.t.disc, "type.sequence"=type.sequence,
              "dat"=dat, "stats"=stats, "type.def"=type.def))
}
