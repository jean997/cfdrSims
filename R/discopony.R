

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
#'@export
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
  chrnum <- as.numeric(strsplit_helper(df$chr, "chr", 2))
  df$chunk <- strsplit_helper(df$file, ".txt", 1, fixed=TRUE )
  df$chunk <- strsplit_helper(df$chunk, "chunk", 2, fixed=TRUE)
  o <- order(chrnum, as.numeric(df$chunk))
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
      fn <- paste0("discopony_output/", c, "/", c, "_maxes.upd.RData")
      ff <- getobj(fn)
      my.names <- sapply(ff, FUN=function(f){f$file})
    }
    k <- which(my.names==df$file[i])
    iv <- name_clusters_merged(x=ff[[k]]$ys, z=df$z[i], z0=0.3*0.9)
    nn <- nrow(iv)
    stopifnot(nn==df$R[i])
    iv <- matrix(ff[[k]]$pos[as.matrix(iv)], ncol=2, byrow=FALSE)
    iv[,2] <- iv[,2]+1
    ivls[j:(j+nn-1), c(3, 4)] <- iv
    j <- j + nn
  }

  ivls <- data.frame(ivls, stringsAsFactors=FALSE)
  names(ivls) <- c("chr", "chunk", "start", "stop", "z")
  for(j in 2:5) ivls[,j] <- as.numeric(ivls[,j])
  ivls$length <- ivls$stop-ivls$start
  return(list("df"=df, "ivls"=ivls))
}


dnase1_test_windows <- function(dat.file, pheno.file, maxit=50){
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
    return(c(beta1, s1, beta1/s1, tp(beta1/s1, df=sum(n)-2)))
  }
  huber_reg <- function(y, labs){
    f <- rlm(y~labs, psi=psi.huber, k=1.345, scale.est="Huber", maxit=50)
    b1 <- summary(f)$coefficients[2, 1]
    s <- summary(f)$coefficients[2, 2]
    if(is.na(s)) return(c(0, 1))
    return(c(b1, s, b1/s, normp(b1/s)))
  }
  tt <- function(y, labs){
    b <- mean(y[labs==1]-mean(y[labs==0]))
    s <- sqrt( var(y[labs==1])/sum(labs==1) + var(y[labs==0])/sum(labs==0))
    return(c(b, s, b/s, tp(b/s, df=length(y)/2-1)))
  }

  dat <- getobj(dat.file)
  X <- read_delim(pheno.file, col_names=FALSE, delim=" ")
  X <- X[match(names(dat)[2:26], X$X1),  ]
  labs <- X$X2


  #Window start stop
  #Four collumns for each stat: beta, se, stat, p-value
  res <- apply(dat, MARGIN=1, FUN=function(x){
    y <- as.numeric(x[2:26])
    c(x[1], x[27], x[28], huber_reg(y, labs), pois_reg(y, labs), tt(y, labs))
  })
  res <- data.frame(t(res))
  nms <- c("Window", "Start", "Stop")
  nms <- c(nms, paste0("Huber", c("Beta", "SE", "Stat", "P")))
  nms <- c(nms, paste0("Pois", c("Beta", "SE", "Stat", "P")))
  nms <- c(nms, paste0("T", c("Beta", "SE", "Stat", "P")))
  names(res) <- nms
  cat("\n")
  return(res)
}

#'@export
dnase1_run_waveqtl <- function(dat.file, pheno.file,
                               window.file, waveQTL_loc, chr, win.range=NULL){
  #Read data. Dat file has first two columns as pos and win
  dat <- read_delim(dat.file, delim=" ")
  X <- read_delim(pheno.file, col_names=FALSE, delim=" ")
  X <- X[match(names(dat)[3:27], X$X1),  ]
  labs <- X$X2

  #Write "genotype" file
  N <- floor(runif(n=1, min=10, max=1e9)) #Random number label
  geno <- c("chr1.1", "A", "G", labs)
  cat(geno, file=paste0("geno_", N, ".txt"))

  #This is a bed file
  win.bound = read_delim(window.file, delim="\t", col_names=FALSE)
  win.bound = win.bound[win.bound$X1==chr,]
  cat(dim(win.bound)[1], " total windows")
  wins = unique(dat$win)
  if(!is.null(win.range)){
    wins=wins[wins >= win.range[1] & wins <= win.range[2]]
  }
  win.bound = as.matrix(win.bound[, 2:3])
  cat(dim(win.bound)[1], " windows to analyze")

  pvals = c()
  for(w in wins){
    cat(w, " ")
    pos = win.bound[w,1]:(win.bound[w,2]-1)
    n = length(pos)
    stopifnot(log2(n)==trunc(log2(n)))
    pdat = dat[dat$win==w & (!dat$pos == win.bound[w,2]), ]
    pheno.dat <- matrix(0, nrow=n, ncol=25)
    ix = match(pdat$pos, pos)
    pheno.dat[ix,] = as.matrix(pdat[, 3:27])
    res <- WaveQTL_preprocess(Data = t(pheno.dat), meanR.thresh = 0)
    f <- tempfile(tmpdir = ".")
    write.table(res$WCs, file=paste0(f, "_pheno.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
    cat(res$filtered.WCs, file=paste0(f, "_use.txt"))
    cmd <- paste0(waveQTL_loc, " -gmode 1 -g geno_", N, ".txt -p ",
                f, "_pheno.txt -u ", f, "_use.txt -o temp", N, " -f ", n, " -numPerm 1000 -fph 2")
    system(cmd)
    pval <- read.table(paste0("output/temp", N, ".fph.pval.txt"), header=TRUE)
    pvals <- c(pvals, pval[3, 1])
    unlink(paste0(f, "_pheno.txt"))
    unlink(paste0(f, "_use.txt"))
    unlink(paste0("output/temp", N, "*"))
  }
  cat("\n")
  win.bound <- win.bound[wins, ]
  ret <- data.frame(cbind(win.bound, pvals))
  names(ret) = c("start", "stop", "pval")
  ret$chr = chr
  return(ret)
}

#'@export
expand_windows <- function(win.file, out.file){
  wins = read_delim(win.file, delim="\t", col_names=FALSE)
  wins$dist = wins$X3-wins$X2
  wins$newDist = 2^ceiling(log2(wins$dist))
  wins$diff = wins$newDist -wins$dist
  wins$start = wins$X2-floor(wins$diff/2)
  wins$stop = wins$X3 + ceiling(wins$diff/2)
  wins.bam = wins[, c("X1", "start", "stop")]
  write.table(wins.bam, file=out.file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

#'@export
strsplit_helper <- function(list, split, field, fixed=FALSE){
  x <- unlist(lapply(list, FUN=function(x){
    unlist(strsplit(x, split, fixed=fixed))[field]}))
}
