

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
#'@param out.file Output file (optional)
#'@param max1.file File with data from running with n.perm=0
#'@param seed Set a seed for running permutations (Required)
#'@param n.perm Number of permutations
#'@param z0 Reference level for merging
#'@param bandwidth Smoother bandwidth.
#'@param maxit Maximum iterations for rlm.
#' @return Saves an RData file with maxes list. Returns the length of that list.
#'@export
discopony_maxes1 <- function(dat.file, pheno.file, s0, zmin,
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

#This function assumes a v. specific file name structure
discopony_pull_regions <- function(results.file, thresh){
  Z <- getobj(results.file)
  if(all(Z$fdr > thresh)){
    cat("No intervals achieve significance")
    return(0)
  }
  ix <- max(which(Z$fdr <= thresh))
  R <- as.numeric(Z$Robs[ix,-1])
  z <- as.numeric(Z$z[ix, -1])
  file.names = names(Z$Robs)[-1]

  df <- data.frame(cbind(z[R > 0],  R[R > 0]), stringsAsFactors = FALSE)
  names(df) <- c("z", "R")
  df$file <- file.names[R > 0]

  ivls <- matrix(nrow=sum(df$R), ncol=5)
  #Chr chunk start stop z
  cat(sum(df$R), " intervals in ", nrow(df), " chunks", "\n")

  df$chr <- strsplit_helper(df$file, "/", 1, fixed=TRUE)
  df$chunk <- strsplit_helper(df$file, ".txt", 1, fixed=TRUE )
  df$chunk <- strsplit_helper(df$chunk, "chunk", 2, fixed=TRUE)
  o <- order(df$chr, as.numeric(df$chunk))
  df <- df[o, ]

  #One row per interval
  ivls[,1] <- rep(df$chr, df$R)
  ivls[,2] <- rep(df$chunk, df$R)
  ivls[,5] <- rep(df$z, df$R)
  j <- 1
  c <- ""
  for(i in 1:nrow(df)){
    if(df$chr[i]!=c){
      c <- df$chr[i]
      cat(c, "..")
      fn <- paste0(c, "/", c, "_maxes.RData")
      ff <- getobj(fn)
      my.names <- sapply(ff, FUN=function(f){f$file})
    }
    k <- which(my.names==df$file[i])
    iv <- name_clusters_merged(x=ff[[k]]$ys, z=df$z[i], z0=0.3*0.9)
    nn <- nrow(iv)
    stopifnot(nn==df$R[i])
    ivls[j:(j+nn-1), c(3, 4)] <- as.matrix(iv)
    j <- j + nn
  }

  ivls <- data.frame(ivls, stringsAsFactors=FALSE)
  names(ivls) <- c("chr", "chunk", "start", "stop", "z")
  for(j in 2:5) ivls[,j] <- as.numeric(ivls[,j])
  ivls$length <- ivls$stop-ivls$start + 1
  return(list("df"=df, "ivls"=ivls))
}


dnase1_test_windows <- function(dat.file, pheno.file, maxit=50, win.range=NULL){
  normp <- function(x){
    if(x < 0) return(2*pnorm(x))
    return(2*pnorm(x, lower.tail=FALSE))
  }
  tp <- function(x, df){
    if(x < 0) return(2*pt(x, df=df))
    return(2*pt(x, df=df, lower.tail=FALSE))
  }
  pois_reg <- function(y, labs, zero.val=1e-11){
    cts <- as.vector( by(y, labs, FUN=function(z){max(zero.val, sum(z))}))
    n <- as.vector( by(labs, labs, FUN=length))
    #poisson
    beta1 <- log(cts[2]/n[2])-log(cts[1]/n[1])
    mu = rep(cts[1]/n[1], sum(n))
    mu[labs==1] = cts[2]/n[2]
    phi = 1/(sum(n) - 2)* sum( (y-mu)^2/mu)
    s1 = sqrt(phi)*sqrt(sum(1/cts))
    return(c(beta1/s1, tp(beta1/s1, df=sum(n)-2)))
  }
  huber_reg <- function(y, labs){
    f <- rlm(y~labs, psi=psi.huber, k=1.345, scale.est="Huber", maxit=50)
    b1 <- summary(f)$coefficients[2, 1]
    s <- summary(f)$coefficients[2, 2]
    return(c(b1/s, normp(b1/s)))
  }
  tt <- function(y, labs){
    f <- t.test(y~labs)
    return(as.vector(c(f$statistic, f$p.value)))
  }
  dat <- read_delim(dat.file, delim=" ")
  X <- read_delim(pheno.file, col_names=FALSE, delim=" ")
  X <- X[match(names(dat)[c(-1, -2)], X[,1]),  ]
  labs <- X[,2]
  wins <- unique(dat$win)
  wins <- as.numeric(wins)
  if(!is.null(win.range)){
    wins <- wins[wins >= win.range[1] & wins <= win.range[2]]
    dat <- dat[dat$win %in% wins,]
  }
  res <- matrix(nrow=length(unique(dat$win)), ncol=9)
  #Window start stop HuberStat HuberP PoisStat PoisP Tstat TP
  for(i in 1:length(wins)){
    cat(i, " ")
    res[i, 1] <- w <- wins[i]
    pos <- dat$pos[dat$win==w]
    res[i, 2] <- min(pos)
    res[i, 3] <- max(pos)
    y <- as.numeric(colSums(dat[dat$win==w, c(-1, -2)]))

    res[i, 4:5] <- huber_reg(y, labs)
    res[i, 6:7] <- pois_reg(y, labs)
    res[i, 8:9] <- tt(y, labs)
  }
  res <- data.frame(res)
  names(res) <- c("Window", "Start", "Stop", "HuberStat", "HuberP", "PoisStat", "PoisP", "TStat", "TP")
  cat("\n")
  return(res)
}





strsplit_helper <- function(list, split, field, fixed=FALSE){
  x <- unlist(lapply(list, FUN=function(x){
    unlist(strsplit(x, split, fixed=fixed))[field]}))
}
