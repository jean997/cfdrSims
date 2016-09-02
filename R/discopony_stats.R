
#' Get stats and permuted stats from dnase 1 data
#'@description Finds peak heights and permuatation peak heights using dnase 1 data.
#'@param dat.file Data file. Should have first
#'column as position and a header giving sample names. One row per base pair.
#'@param pheno.file Phenotype file. First column should correspond to the header of dat.file.
#'@param s0 Variance inflation constant.
#'@param zmin Lower bound on significance thresholds.
#'@param z0 Reference level for merging
#'@param strt Region start base pair (may not be the same as where dat.file starts)
#'@param stp Region stop base pair.
#'@param out.file Output file (optional)
#'@param max1.file File with data from running with n.perm=0
#'@param seed Set a seed for running permutations (Required)
#'@param n.perm Number of permutations
#'@param z0 Reference level for merging
#'@param bandwidth Smoother bandwidth.
#'@param maxit Maximum iterations for rlm.
#' @return Saves an RData file with maxes list. Returns the length of that list.
#'@export
discopony_stats <- function(dat.file, pheno.file, s0, zmin,
                             seed, n.perm, strt=NULL, stp=NULL,
                             out.file=NULL, max1.file=NULL,
                             z0=zmin*0.3, bandwidth=50, maxit=50){

  X <- read_delim(pheno.file, col_names=FALSE, delim=" ")
  n <- nrow(X)

  #Permuted phenotypes
  if(n.perm > 0){
    set.seed(seed)
    perms <- replicate(n=n.perm, expr = {
      sample( X[,2], size=n, replace=FALSE)
    })
  }

  #Read data
  name.root <- unlist(strsplit(dat.file, ".txt"))[1]
  dat <- read_delim(dat.file, delim=" ")
  if(nrow(dat)==0) return(0)
  pos <- dat[,1]
  #Make sure the phenotype is sorted correctly
  X <- X[match(names(dat)[-1], X[,1]),  ]
  #Chunks will have a little extra data to get the smoothing right
  if(!is.null(strt)){
    ix1 <- min(which(pos >=strt))
  }else{
    ix1 <- 1
    strt <- pos[1]
  }
  if(!is.null(stp)){
    ix2 <- max(which(pos <= stp))
  }else{
    ix2 <- length(pos)
    stp <- pos[ix2]
  }
  len <- stp-strt + 1

  #Calculate statistics
  if(is.null(max1.file)){
    y <- huber_stats2(Y=dat[, -1], labs=X[,2],s0=s0, maxit=maxit)
    ys <- ksmooth_0(x=pos, y=y, bandwidth = bandwidth)[ix1:ix2]
  }else{
    m1 <- getobj(max1.file)
    ys <- m1$ys
  }
  if(all(abs(ys) < zmin)){
    cat("No clusters exceed ", zmin, "\n")
    #unlink(file.name)
    return(0)
  }

  #Find peak heights
  #Intervals defined by z0:
  q0 <-rle( abs(ys) > z0 )
  p0 <- length(q0$lengths)
  ivls <- cbind(c(1, cumsum(q0$lengths)[-p0]+1)[q0$values], (cumsum(q0$lengths))[q0$values])
  max1 <- apply(ivls, MARGIN=1, FUN=function(iv){ max(abs(ys)[iv[1]:iv[2]])})
  max1 <- max1[max1 >= zmin]
  if(n.perm==0){
    R <- list("max1"=max1,"file"=dat.file, "nbp"=len, "ys"=ys, "pos"=pos, "z0"=z0, "zmin"=zmin)
    return(R)
  }

  #Permutation test statistics and peak heights
  cat("Calculating permutation peak heights.\n")
  max.perm <- apply(perms, MARGIN=2, FUN=function(l){
    yy <- huber_stats2(dat[,-1], labs=l, s0=s0, maxit=maxit)
    yys <- ksmooth_0(x=pos, y=yy, bandwidth = bandwidth)[ix1:ix2]
    if(all(abs(yys) <= z0)) return(c())
    q0 <-rle( abs(yys) > z0 )
    p0 <- length(q0$lengths)
    ivls <- cbind(c(1, cumsum(q0$lengths)[-p0]+1)[q0$values], (cumsum(q0$lengths))[q0$values])
    apply(ivls, MARGIN=1, FUN=function(iv){ max(abs(yys)[iv[1]:iv[2]])})
  })
  m <- sort(unlist(max.perm), decreasing=TRUE)

  mx <- cbind(m, (1:length(m))/(n.perm*len))

  if(all(m < zmin)){
    mx <- cbind(zmin, 0)
  }else if(sum(m >=zmin) < 5){
    mx <- mx[1:5, ]
  }else{
    mx <- mx[m >= zmin,]
  }

  if(is.null(out.file)){
    out.file <- paste0(name.root, "_mx.RData")
  }

  R <- list("max1"=max1, "mx"=mx, "file"=dat.file,
            "z0"=z0, "zmin"=zmin,
            "nbp"=len, "ys"=ys, "pos"=pos)
  save(R, file=out.file)
  return(nrow(mx))
}
