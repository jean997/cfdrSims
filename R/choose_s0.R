

choose_s0_z0 <- function(Y, labs, pos, maxit=50, z0_quantile=0.9, bandwidth=20){

  x <- cfdrSims:::huber_helper(Y=Y, labs=labs, maxit=maxit)
  salpha <- c(0, as.numeric(quantile(x[2,], probs=seq(0, 1, by=0.05))))
  nn <- length(salpha)
  ix <- sapply(x[2,], FUN=function(w){ sum(w >= salpha[-nn])})
  cvs <- c()
  v <- as.numeric(by(data=x[1,], INDICES = ix, FUN = mad))
  cvs <- c(sd(v)/mean(v))
  tol <- min(x[2,][x[2,] > 0])/10

  fct <- function(s0, Y, labs, ix,  maxit){
    cat(s0, " ")
    xx <- cfdrSims:::huber_stats2(Y, labs, s0=s0, maxit=maxit)
    v <- as.numeric(by(data=xx, INDICES = ix, FUN = mad))
    return(sd(v)/mean(v))
  }

  W <- optimize(fct, interval=range(x[2,]), Y=Y,
                labs=labs, ix=ix, maxit=maxit, tol=tol)
  xx <- cfdrSims:::huber_stats2(Y, labs, s0=W$minimum, maxit=maxit)
  xs <- ksmooth(x=pos, y=xx, bandwidth=bandwidth)$y
  z0 <- as.numeric(quantile(xs, probs=z0_quantile))
  return(list("s0"=W$minimum, "z0"=z0))
}

dnase_choose_s0_z0 <- function(file.name, pheno.file, seed, has.stats=FALSE, n.select=1e4){
  set.seed(seed)
  dat <- read_table(file.name)
  n <- nrow(dat)
  p <- ncol(dat)
  if(has.stats){
    dat <- dat[,-p]
    p <- p-1
  }
  ix <- sample(1:n, size=n.select, replace=FALSE)
  dat <- dat[ix,]
  X <- read_table(pheno.file, col_names=FALSE)
  x <- X[match(names(dat)[-1], X[,1]),2]
  s= choose_s0_z0(Y=dat, labs = x)
  return(s)
}
