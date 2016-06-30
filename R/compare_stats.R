
compare_stats <- function(group1.func, group2.func, n.rep=1000, sample.size=c(50, 50)){
  #Methods to comare are
  #lm, poisson, quasipoisson, negative binomial,
  #huber, poisson-huber
  #nms <- c("lm", "pois", "qpois", "nb",
   #        "huber-lm","huber-pois", "huber2", "huber2-6.5")
  nms <- c("lm", "pois", "qpois", "nb", "huber")

  ps <- matrix(nrow=n.rep, ncol=length(nms))
  betas <- matrix(nrow=n.rep, ncol=length(nms))

  labs <- rep(c(0, 1), sample.size)
  Y <- matrix(nrow=n.rep, ncol=sum(sample.size))
  for(j in 1:n.rep){
    Y[j, labs==0] <- group1.func(n=sample.size[1])
    Y[j, labs==1] <- group2.func(n=sample.size[2])

    p.values <- c()


    #lm
    f <- lm(Y[j, ]~ labs)
    betas[j, 1] <- f$coefficients[2]
    ps[j, 1] <- summary(f)$coefficients[2, 4]
    #Poisson regression: pois
    f <- glm(Y[j, ]~ labs, family=poisson(link="log"))
    betas[j, 2] <- f$coefficients[2]
    ps[j, 2] <- summary(f)$coefficients[2, 4]
    #Quasi-Poisson
    f <- glm(Y[j, ]~ labs, family=quasipoisson(link="log"))
    betas[j, 3] <- f$coefficients[2]
    ps[j, 3] <- summary(f)$coefficients[2, 4]
    #Negative Binomial
    f <- glm.nb(Y[j, ]~ labs)
    betas[j, 4] <- f$coefficients[2]
    ps[j, 4] <- summary(f)$coefficients[2, 4]


    #Huber
    #f <- lmrob(Y[j,]~labs)
    #betas[j, 5] <- f$coefficients[2]
    #ps[j, 5] <- summary(f)$coefficients[2, 4]

    #Huber - poisson
    #f <- glmrob(Y[j, ]~ labs, family=poisson(link="log"))
    #betas[j, 6] <- f$coefficients[2]
    #ps[j, 6] <- summary(f)$coefficients[2, 4]

    #Huber -2
    f <- rlm(Y[j,]~labs)
    betas[j, 5] <- f$coefficients[2]
    #ps[j, 7] <- get.pt(summary(f)$coefficients[2,3], df=sum(sample.size)-2)
    ps[j, 5] <- get.p(summary(f)$coefficients[2,3])
    #v <- huber_var_eq6.5(f)
    #ps[j, 8] <- get.pt(betas[j, 7]/sqrt(v[2, 2]), df=sum(sample.size)-2)
    #ps[j, 8] <- get.p(betas[j, 7]/sqrt(v[2, 2]))

  }
  ps <- data.frame(ps)
  names(ps) <- nms
  betas <- data.frame(betas)
  names(betas)<- nms[-8]
  return(list("p"=ps, "beta"=betas, "Y"=Y))
}

#Seed 209601188
run_compare <- function(seed=209601188, alphas=c(0.025, 0.05, 0.1), sample.size=c(50, 50), n.rep=500){
  set.seed(seed)
  t1err.tab <- list()
  power.tab <- list()
  ZE <- list()
  ZD <-list()
  log.methods <- c("qpois", "nb")
  id.methods <- c("lm", "huber")
  methods <- c(log.methods, id.methods)

  settings <- c()
  for(i in 1:length(alphas)){
    t1err.tab[[i]] <- power.tab[[i]] <- matrix(nrow=5, ncol=length(methods))
  }

  #Setting 1 - Poisson
  i <- 1
  cat(i, "\n")
  g1 <- function(n){rpois(n=n, lambda=3)}
  g2 <- function(n){rpois(n=n, lambda=4)}
  ZD[[i]] <- cfdrSims:::compare_stats(g1, g2, n.rep, sample.size)
  ZE[[i]] <- cfdrSims:::compare_stats(g1, g1, n.rep, sample.size)
  settings  <- c(settings, "Poisson")

  #Setting 2 - Poisson with outliers
  i <- 2
  cat(i, "\n")
  g1 <- function(n){
    p <- sample(0:1, prob=c(0.95, 0.05), size=n, replace=TRUE)
    lams <- rep(3, n)
    lams[p==1] <- 15
    rpois(n=n, lambda=lams)
  }
  g2 <- function(n){
    p <- sample(0:1, prob=c(0.95, 0.05), size=n, replace=TRUE)
    lams <- rep(4, n)
    lams[p==1] <- 15
    rpois(n=n, lambda=lams)
  }
  ZD[[i]] <- compare_stats(g1, g2, n.rep, sample.size)
  ZE[[i]] <- compare_stats(g1, g1, n.rep, sample.size)
  settings <- c(settings, "Poisson-Outliers")

  #Setting 5 Multinomial - Poisson - signal in large rates
  i <- 3
  cat(i, "\n")
  g1 <- function(n){
    lambda <- sample(c(3, 15), prob =c(0.95, 0.05), size=n, replace=TRUE)
    rpois(n=n, lambda=lambda)}
  g2 <- function(n){
    lambda <- sample(c(3,  15), prob =c(0.9, 0.1), size=n, replace=TRUE)
    rpois(n=n, lambda=lambda)}
  ZD[[i]] <- compare_stats(g1, g2, n.rep, sample.size)
  ZE[[i]] <- compare_stats(g1, g1, n.rep, sample.size)
  settings <- c(settings, "Poisson-Outliers2")



  #Setting 3 - Exponential Poisson
  i <- 4
  cat(i, "\n")
  g1 <- function(n){lambda = rexp(n, rate=1/3); rpois(n=n, lambda=lambda)}
  g2 <- function(n){lambda = rexp(n, rate=1/4); rpois(n=n, lambda=lambda)}
  ZD[[i]] <- compare_stats(g1, g2, n.rep, sample.size)
  ZE[[i]] <- compare_stats(g1, g1, n.rep, sample.size)
  settings <- c(settings, "Exponential-Poisson")

  #Setting 4 Multinomial - Poisson - signal not in large rates
  i <- 5
  cat(i, "\n")
  g1 <- function(n){
    lambda <- sample(c(1.5, 3, 7, 15), prob =c(0.55, 0.3, 0.1, 0.05), size=n, replace=TRUE)
    rpois(n=n, lambda=lambda)}
  g2 <- function(n){
    lambda <- sample(c(1.5, 3, 7, 15), prob =c(0.2, 0.45, 0.25, 0.1), size=n, replace=TRUE)
    rpois(n=n, lambda=lambda)}
  ZD[[i]] <- compare_stats(g1, g2, n.rep, sample.size)
  ZE[[i]] <- compare_stats(g1, g1, n.rep, sample.size)
  settings <- c(settings, "Multinomial-Poisson")


  for(i in 1:length(settings)){
    for(m in 1:length(methods)){
      j <- which(names(ZD[[i]]$p)==methods[m])
      for(k in 1:length(alphas)){
        t1err.tab[[k]][i, m] <- mean(ZE[[i]]$p[,j] < alphas[k], na.rm=TRUE)
        power.tab[[k]][i, m] <- mean(ZD[[i]]$p[,j] < alphas[k], na.rm=TRUE)
      }
    }
  }

  for(k in 1:length(alphas)){
    t1err.tab[[k]] <- data.frame(t1err.tab[[k]])
    power.tab[[k]] <- data.frame(power.tab[[k]])
    names(power.tab[[k]]) <- names(t1err.tab[[k]]) <- methods
    rownames(power.tab[[k]]) <- rownames(t1err.tab[[k]]) <- settings
  }
  return(list("t1E" = t1err.tab, "power"=power.tab, "ZD"=ZD, "ZE"=ZE))
}


make_roc_curves <- function(ZD, ZE, which.setting, which.methods=NULL, cols=NULL, main="", ltys=NULL){
  d <- ZD[[which.setting]]$p
  e <- ZE[[which.setting]]$p
  if(is.null(which.methods)) which.methods <- 1:ncol(d)
  if(is.null(ltys)) ltys <- rep(1, length(which.methods))
  n <- nrow(d)
  if(is.null(cols)) cols <- 1:length(which.methods)
  for(j in 1:length(which.methods)){
    i <- which.methods[j]
    pred <- prediction(-log(c(d[,i], e[,i])), labels = rep(c(1,0), each=n))
    perf <- performance(pred, "tpr", "fpr")
    if(j==1) add=FALSE
      else add=TRUE
    plot(perf, add=add, col=cols[j], main=main, lty=ltys[j], lwd=1.5)
  }
  legend("bottomright", legend=names(d)[which.methods], col=cols, lty=ltys, cex=2)
}



get.res <- function(m1, m2, Z){
  log.methods <- c("pois", "qpois", "nb", "huber-pois")
  id.methods <- c("lm", "huber-lm", "huber2")
  k <- nrow(Z$beta)
  res <- matrix(nrow=7, ncol=5)
  epois <- var(Z$beta[, names(Z$beta)=="pois"])
  elm <- var(Z$beta[, names(Z$beta)=="lm"])
  ct <- 1
  for(m in log.methods){
    i <- which(names(Z$beta)==m)
    bias <- mean(Z$beta[,i] - (log(m2/m1)))
    rel.eff <- epois/var(Z$beta[,i])
    p1 <- sum(Z$p[,i] < 0.025, na.rm=TRUE)/k
    p2 <- sum(Z$p[,i] < 0.05, na.rm=TRUE)/k
    p3 <- sum(Z$p[,i] < 0.1, na.rm=TRUE)/k
    res[ct, ] <- c(bias, rel.eff, p1, p2, p3)
    ct <- ct + 1
  }
  for(m in id.methods){
    i <- which(names(Z$beta)==m)
    bias <- mean(Z$beta[,i] - (m2-m1))
    rel.eff <- elm/var(Z$beta[,i])
    p1 <- sum(Z$p[,i] < 0.025, na.rm=TRUE)/k
    p2 <- sum(Z$p[,i] < 0.05, na.rm=TRUE)/k
    p3 <- sum(Z$p[,i] < 0.1, na.rm=TRUE)/k
    res[ct, ] <- c(bias, rel.eff, p1, p2, p3)
    ct <- ct + 1
  }
  res <- data.frame(res)
  names(res) <- c("bias", "relEff", "p<0.025", "p<0.05", "p<0.1")
  rownames(res) <- c(log.methods, id.methods)
  return(res)
}


get.p <- function(x){
  if(x < 0) return(2*pnorm(x))
  return(2*pnorm(x, lower.tail=FALSE))
}

get.pt <- function(x, df){
  if(x < 0) return(2*pt(x, df=df))
  return(2*pt(x, df=df, lower.tail=FALSE))
}
