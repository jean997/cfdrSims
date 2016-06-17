

#' Get peak heights and permuted peak heights from dnase 1 data
#'@description Finds peak heights and permuatation peak heights using dnase 1 data.
#'@param dat.file Data file. Should have first
#'column as position and a header giving sample names. One row per base pair.
#'@param pheno.file Phenotype file. First column should correspond to the header of dat.file.
#'@param s0 Variance inflation constant.
#'@param zmin Lower bound on significance thresholds.
#'@param z0 Reference level for merging
#'@param strt Region start base pair (may not be the same as where dat.file starts)
#'@param stp Region stop base pair.
#'@param seed Set a seed for running permutations (Required)
#'@param n.perm Number of permutations
#'@param z0 Reference level for merging
#'@param bandwidth Smoother bandwidth.
#'@param maxit Maximum iterations for rlm.
#' @return Saves an RData file with maxes list. Returns the length of that list.
#'@export
discopony_maxes1 <- function(dat.file, pheno.file, s0, zmin,
                      seed, n.perm, strt=NULL, stp=NULL,
                      z0=zmin*0.3, bandwidth=50, maxit=50){

  X <- read_delim(pheno.file, col_names=FALSE, delim=" ")
  n <- nrow(X)

  #Permuted phenotypes
  set.seed(seed)
  perms <- replicate(n=n.perm, expr = {
    sample( X[,2], size=n, replace=FALSE)
  })

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
  #Calculate statistics
  y <- huber_stats2(Y=dat[, -1], labs=X[,2],s0=s0, maxit=maxit)
  ys <- ksmooth_0(x=pos, y=y, bandwidth = bandwidth)[ix1:ix2]
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
  len <- stp-strt + 1
  mx <- cbind(m, (1:length(m))/(n.perm*len))

  if(all(m < zmin)){
    mx <- cbind(zmin, 0)
  }else if(sum(m >=zmin) < 5){
    mx <- mx[1:5, ]
  }else{
    mx <- mx[m >= zmin,]
  }
  file.name <- paste0(name.root, "_mx.RData")

  R <- list("max1"=max1, "mx"=mx, "file"=dat.file, "nbp"=len, "ys"=ys, "pos"=pos)
  save(R, file=file.name)
  return(nrow(mx))
}




#' Find thresholds for a range of lambda values. Calculate FDR.
#'@description Find thresholds for a range of lambda values. Calculate FDR.
#'@param file.list List of files. Each file should contain a list of objects produced by
#'discopony_maxes1
#'@param zmin Lower bound on significance thresholds.
#'@param nlam Number of lambda values to consider
#' @return A list with items z, Robs, and fdr.
#'@export
discopony_choose_z <- function(file.list, zmin, nlam){
  log.lambda.min <- Inf
  log.lambda.max <- -Inf
  nbp <- 0
  n.chunk <- 0
  names <- c()
  for(f in file.list){
    R <- getobj(f)
    for(i in 1:length(R)){
      log.lambda.min <- min(log.lambda.min, log10(R[[i]]$mx[,2][R[[i]]$mx[,2] > 0 & R[[i]]$mx[,1] >=zmin]))
      log.lambda.max <- max(log.lambda.max, log10(R[[i]]$mx[,2][R[[i]]$mx[,1] >=zmin]))
      names <- c(names, R[[i]]$file)
      nbp <- nbp + R[[i]]$nbp
    }
    n.chunk <- n.chunk + length(R)
  }
  lams <- seq(log.lambda.min, log.lambda.max, length.out=nlam)

  n.seg <- length(file.list)
  z <- matrix(nrow=nlam, ncol=n.chunk+1)
  Robs <- matrix(nrow=nlam, ncol=n.chunk + 1)
  z[,1] <- Robs[, 1] <- lams
  ct <- 1
  for(f in file.list){
    R <- getobj(f)
    for(i in 1:length(R)){
      if(nrow(R[[i]]$mx)==1 & R[[i]]$mx[1, 2]==0){
        z[, ct+1] = zmin
        Robs[, ct + 1] = sum(R[[i]]$max1 >= zmin)
      }else{
        z[, ct + 1] <- approx(y=R[[i]]$mx[,1], x=log10(R[[i]]$mx[,2]),
                              xout=lams, yright=zmin, yleft=Inf)$y
        Robs[, ct + 1] <- sapply(z[, ct+1], FUN=function(zz){sum(R[[i]]$max1 > zz)})
      }
      ct = ct + 1
    }
  }
  z <- data.frame(z)
  names(z)<- c("lambda", names)
  Robs <- data.frame(Robs)
  names(Robs) <- c("lambda", names)
  fdr <- (10^(Robs[,1]))*nbp/rowSums(Robs[, -1, drop=FALSE])
  ret <- list("Robs"=Robs, "z"=z, "fdr"=fdr)
  return(ret)
}

