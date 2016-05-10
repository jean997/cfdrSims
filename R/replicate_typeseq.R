
#' Replicate one type sequence many times
#'@description Generate data from a type sequence to estimate z
#'@param type.sequence Vector of integers in 1:6 giving the type of each subregion.
#'@param sample.size Vector of length 2 giving sample sizes
#'@param seed Optional - set a seed.
#'@param n.rep Number of simulations
#'@param s0 Vector of length 3 giving additional variance to add to Poisson, Huber and t-test statistics.
#'@param n.perms Number of permutations to use to estimate lambda.
#' @return A list of statistic matrices, one for each test type
#'@export
replicate_typeseq <- function(type.sequence, file=NULL, sample.size=c(20, 20),
                              seed=NULL, n.rep=500, s0=c(0, 0, 0),
                               type.def=NULL,
                              stat.names = c("Poisson", "Huber", "t-test")){

  if(!is.null(seed)) set.seed(seed)
  if(!is.null(type.def)) stopifnot(names(type.def) ==c("p1", "p2"))
  labs <- rep(c(0, 1), sample.size)


  D <- replicate(n=n.rep,
                 expr=sample_data(type.sequence, sample.size, type.def)$dat)

  stats <- list()

  for(i in 1:length(stat.names)){
    if(stat.names[i]=="Poisson") stats[[i]] <- apply(D, MARGIN=c(3), FUN=function(x){
      pois_regression(Y=x, labs=labs, s0=s0[i])})
    else if(stat.names[i]=="Huber") stats[[i]]<- apply(D, MARGIN=c(3), FUN=function(x){
      huber_stats(Y=x, labs=labs, s0=s0[i])})
    else if(stat.names[i]=="t-test") stats[[i]] <- apply(D, MARGIN=c(3), FUN=function(x){
      t_stats(dat=x, labs=labs, s0=s0[i])})
  }

  names(stats) <- stat.names


  r <- length(type.sequence)
  q <- sum(type.sequence %in% 4:6)
  p <- 200 * length(type.sequence)
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
  }else{
    signal <- Intervals()
  }

  R <- list("stats"=stats, "type.sequence"=type.sequence, "signal"=signal)
  if(!is.null(file)) save(R, file=file)
  return(R)
}
