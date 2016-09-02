lilexample_profile <- function(){
  yy <- cfdrSims:::gen_profile(5)
  yy <- c(yy,  cfdrSims:::gen_profile(3))
  yy <- c(yy,  cfdrSims:::gen_profile(7))
  #yy <- c(yy,  cfdrSims:::gen_profile(c(3.2, 4.8)))
  yy <- c(yy,  cfdrSims:::gen_profile(4.8))
  df <- data.frame(pos=1:800, y=yy)

  p <- 800; k <- 4
  w_64b <- cbind( 69 + 200*(1:k -1), 132 + 200*(1:k -1))
  informed <- data.frame(w_64b)
  names(informed) = c("x", "xend")
  informed$type <- "informed"
  informed$y <- 1
  w_64e <- cbind(seq(1, p - (p%%64), by=64), seq(64, p, by=64))
  naive <- data.frame(w_64e)
  names(naive)=c("x", "xend")
  naive$type="naive"
  naive$y <-rep(c(0.5, 0.6), 6)
  allbins <- rbind(informed, naive)

  plot <- ggplot(df) + geom_line(aes(x=pos, y=yy)) +
    geom_segment(data=allbins, aes(x=x, xend=xend, y=y, yend=y, col=type), lwd=1.5) +
    scale_color_manual(labels=c("Informed Bins", "Naive Bins"), values=c("blue", "darkorange2")) +
    scale_y_continuous(breaks=c(1.5, 2:7)) + ylab(expression("Profile ("~gamma~")")) +
    theme_bw(18) + theme(axis.title.x=element_blank(),
                         panel.grid=element_blank(),
                         legend.title=element_blank(),
                         legend.position=c(0.1, 0.85),
                         legend.key = element_blank(),
                         legend.text = element_text(size=15))
  ggsave(plot, file="~/Dropbox/Thesis/img/example_profile2.png", height=4.5, width=15, units="in", dpi=300)


  plot_talk <- ggplot(df) + geom_line(aes(x=pos, y=yy)) +
    geom_segment(data=allbins, aes(x=x, xend=xend, y=y, yend=y, col=type), lwd=1.5) +
    scale_color_manual(labels=c("Informed Bins", "Naive Bins"), values=c("blue", "darkorange2")) +
    scale_y_continuous(breaks=c(1.5, 2:7)) + ylab(expression("Profile ("~gamma~")")) +
    xlab("Position") +
    theme_bw(18) + theme(#axis.title.x=element_blank(),
                         panel.grid=element_blank(),
                         legend.title=element_blank(),
                         legend.position=c(0.2, 0.85),
                         legend.key = element_blank(),
                         legend.key.width=unit(0.3, "in"),
                         legend.background = element_blank(),
                         legend.text = element_text(size=15))
  ggsave(plot_talk, file="~/Dropbox/Thesis/final_talk/img/example_profile2.png",
         height=3, width=7, units="in", dpi=300)

}


plot_lambda <- function(){
  R <- getobj("cont1_11_fret.RData")

  dfz1 <- data.frame(R$cl[[2]][[1]]$z)
  dfz1$R <- R$cl[[2]][[1]]$R[,2]
  names(dfz1) = c("lambda", "z", "r")
  dfz1$lambda <- 12000 * 10^dfz1$lambda
  dfz1$loverR <- dfz1$lambda/dfz1$r
  dfz2 <- data.frame(R$cl[[2]][[2]]$z)
  names(dfz2) <- c("lambda", "z1", "z2")
  dfz2$lambda <- 12000 * 10^dfz2$lambda
  dfz2_long <- gather(dfz2, "interval", "z", -lambda)

  labs <- c( expression(widehat(lambda)[1]~"(z)"),
             expression(widehat(lambda)[2]~"(z)"))


  h1 <- ggplot(dfz1) + geom_line(aes(x=z, y=lambda)) +
  xlab("z") + ylab(expression(lambda)) +
  scale_y_log10()+
  #annotate("text", x=2.5, y=1, label="widehat(lambda)~'(z)'", parse=TRUE, size=7)+
  theme_bw(18) + theme( panel.grid=element_blank())
  ggsave(h1, file="~/Dropbox/Thesis/img/lambda_z2_1.png", height=4.5, width=4.5, units="in", dpi=300)

  h00_talk <- ggplot(dfz1) + geom_line(aes(x=z, y=lambda)) +
    xlab(expression("z")) + ylab(expression(lambda~"(z)"=="E["~V[z]~"]")) +
    scale_y_log10() +
    theme_bw(18) + theme( panel.grid=element_blank(), panel.border = element_blank(),
                          axis.line.x=element_line(), #axis.text=element_blank(),
                          axis.line.y=element_line())
  ggsave(h00_talk, file="~/Dropbox/Thesis/final_talk/img/lambda_00.png", height=4.5, width=4.5, units="in", dpi=300)


  #Plot lambda/R v z
  ix <- which(dfz1$r> 0)
  alph=0.1
  ixx <- max(which(dfz1$loverR <= alph))
  ztilde <- approx(x=dfz1$loverR[(ixx):(ixx+1)], y=dfz1$z[(ixx):(ixx+1)], xout=alph)$y
  h01_talk <- ggplot(dfz1[ix,]) + geom_line(aes(x=z, y=loverR)) +
    geom_point(aes(x=ztilde, y=alph), size=4) +
    geom_segment(aes(x=ztilde, xend=ztilde, y=-Inf, yend=alph), lty=2) +
    xlab(expression("z")) + ylab(expression(lambda~"(z)/"~R[z])) +
    geom_hline(yintercept = 0.1, lty=2) +
    #geom_vline(xintercept = zstar, lty=3) + geom_hline(yintercept = lstar, lty=3) +
    scale_y_continuous(breaks=c(0, 0.1, 0.25, 0.5, 0.75, 1), labels=c(0, expression(alpha==0.1), 0.25, 0.5, 0.75, 0.1))+
    scale_x_continuous(breaks=c(0.5, 1, ztilde, 2), labels=c(0.5, 1, expression(tilde(z)), 2)) +
    #annotate("text", x=2.5, y=1, label="widehat(lambda)~'(z)'", parse=TRUE, size=7)+
    theme_bw(18) + theme( panel.grid=element_blank(), panel.border = element_blank(),
                          axis.line.x=element_line(), #axis.text=element_blank(),
                          axis.line.y=element_line())
  ggsave(h01_talk, file="~/Dropbox/Thesis/final_talk/img/lambda_01.png", height=4.5, width=4.5, units="in", dpi=300)



  L <- approx(x=dfz1$z, y=dfz1$lambda, xout=ztilde)$y
  z0 <- ztilde
  zstar <- 0.9
  lstar <- approx(x=dfz1$z, y=dfz1$lambda, xout=zstar)$y
  line.dat1 = data.frame(x = c(-Inf, z0, z0),  y = c(L, L, 0))

  h1_talk <- ggplot(dfz1) + geom_line(aes(x=z, y=lambda)) +
    geom_line(data=line.dat1, aes(x=x, y=y), lty=2) +
    geom_point(aes(x=z0, y=L), size=4) +
    xlab(expression("z")) + ylab(expression(lambda~"(z)"=="E["~V[z]~"]")) +
    geom_vline(xintercept = zstar, lty=3) + geom_hline(yintercept = lstar, lty=3) +
    scale_y_log10(breaks=c(lstar, L), labels=c(expression(lambda^"*"), expression(Lambda)))+
    scale_x_continuous(breaks=c(zstar, z0), labels=c(expression("z*"), expression(tilde("z")))) +
    #annotate("text", x=2.5, y=1, label="widehat(lambda)~'(z)'", parse=TRUE, size=7)+
    theme_bw(18) + theme( panel.grid=element_blank(), panel.border = element_blank(),
                          axis.line.x=element_line(), #axis.text=element_blank(),
                          axis.line.y=element_line())
  ggsave(h1_talk, file="~/Dropbox/Thesis/final_talk/img/lambda_1.png", height=4.5, width=4.5, units="in", dpi=300)



  h2 <- ggplot(dfz2_long) + geom_line(aes(x=z, y=lambda, group=interval,
                                        color=interval, linetype=interval)) +
  xlab("z") +
  scale_y_log10() + ylab("") +
  theme_bw(18) + theme(panel.grid=element_blank(),
                     legend.title=element_blank(),
                     legend.position=c(0.75, 0.8)) +
  scale_color_manual(labels=labs, values=c("blue", "darkorange2")) +
  scale_linetype_manual(labels=labs, values=c(1, 2))

  ggsave(h2, file="~/Dropbox/Thesis/img/lambda_z2_2.png",
       height=4.5, width=4.5, units="in", dpi=300)


  #labs <- c("Interval 1", "Interval 2")

  L <- L/2
  z1 <- approx(x=dfz2$lambda, y=dfz2$z1, xout=L)$y
  z2 <- approx(x=dfz2$lambda, y=dfz2$z2, xout=L)$y

  line.dat1 = data.frame(x = c(-Inf, z2, z2),  y = c(L, L, 0))
  line.dat2 <- data.frame(x = c(z1, z1),  y = c(0, L))
  point.dat <- data.frame("lam"=rep(L, 2), "z"=c(z1, z2), interval=factor(c("z1", "z2")))
  h2_talk <- ggplot(dfz2_long) + geom_line(aes(x=z, y=lambda, group=interval,
                                               color=interval, linetype=interval), lwd=1.2) +
    geom_line(data=line.dat1, aes(x=x, y=y), lty=2) +
    geom_line(data=line.dat2, aes(x=x, y=y), lty=2) +
    geom_point(data=point.dat, aes(x=z, y=lam, color=interval), size=4) +
    scale_color_manual(labels=labs, values=c("blue", "chartreuse3")) +
    scale_linetype_manual(labels=labs, values=c(1, 1)) +
    xlab(expression("z")) +
    scale_y_log10(breaks=L, labels=expression(Lambda~"/2")) +
    scale_x_continuous(breaks=c(z1, z2),  labels=c(expression(tilde("z")[1]), expression(tilde("z")[2]))) +
    ylab(expression(lambda~"(z)"=="E["~V[z]~"]")) +
    theme_bw(18) + theme(panel.grid=element_blank(),
                         legend.title=element_blank(),
                         panel.border = element_blank(),
                         axis.line.x=element_line(), #axis.text=element_blank(),
                         axis.line.y=element_line(),
                         legend.key = element_blank(),
                         legend.key.width=unit(0.3, "in"),
                         legend.position=c(0.6, 0.8))

  ggsave(h2_talk, file="~/Dropbox/Thesis/final_talk/img/lambda_2.png",
         height=4.5, width=4.5, units="in", dpi=300)

  h20_talk <- ggplot(dfz2_long) + geom_line(aes(x=z, y=lambda, group=interval,
                                               color=interval, linetype=interval), lwd=1.2) +
    geom_line(data=line.dat1, aes(x=x, y=y), lty=2) +
    geom_line(data=line.dat2, aes(x=x, y=y), lty=2) +
    geom_point(data=point.dat, aes(x=z, y=lam, color=interval), size=4) +
    scale_color_manual(labels=labs, values=c("blue", "chartreuse3")) +
    scale_linetype_manual(labels=labs, values=c(1, 1)) +
    xlab(expression("z")) +
    scale_y_log10(breaks=L, labels=expression(frac(lambda["total"], 2))) +
    scale_x_continuous(breaks=c(z1, z2),  labels=c(expression(z[1]), expression(z[2]))) +
    ylab(expression(lambda~"(z)"=="E["~V[z]~"]")) +
    theme_bw(18) + theme(panel.grid=element_blank(),
                         legend.title=element_blank(),
                         panel.border = element_blank(),
                         axis.line.x=element_line(), #axis.text=element_blank(),
                         axis.line.y=element_line(),
                         legend.key = element_blank(),
                         legend.key.width=unit(0.3, "in"),
                         legend.position=c(0.6, 0.8))

  ggsave(h20_talk, file="~/Dropbox/Thesis/final_talk/img/lambda_20.png",
         height=4.5, width=4.5, units="in", dpi=300)





}

plot_lambda_new <- function(){
  mu <-rep(c(0, 1.5, 0, 1.5, 0, 1.5, 0), c(50, 10, 50, 10, 80, 10, 40))
  q <- rle(mu)
  ll <- cumsum(q$lengths)
  ll <- ll[-length(ll)]
  signal <- Intervals(matrix(ll, ncol=2, byrow=TRUE))
  set.seed(1111111)
  n.rep <- 500
  #n.rep <- 10
  Y <- replicate(n=n.rep, expr = {
    y <- rnorm(n=length(mu), mean=mu)
    ksmooth(x=1:length(mu), y=y, x.points = 1:length(mu), bandwidth = 10)$y
  })

  nY <- Y
  for(i in 1:nrow(signal)) nY[signal[i, 1]:signal[i, 2], ] <- 0
  z <- seq(quantile(abs(nY), probs = 0.8), max(abs(nY)), length.out=200)

  lambda <- apply(Y, MARGIN=2,  FUN=function(ys){
    sapply(z, FUN=function(zz){
      if(all(abs(ys) > zz)) return(0)
      q <-rle( abs(ys) > zz )
      p <- length(q$lengths)
      qI <- Intervals(cbind(c(1, cumsum(q$lengths)[-p]+1)[q$values], (cumsum(q$lengths))[q$values] ))
      dd <- distance_to_nearest(qI, signal)
      sum(dd > 0)
    })
  })
  lams <- rowMeans(lambda)
  dfz1 <- data.frame(z=z, lambda=lams)

  h1 <- ggplot(dfz1) + geom_line(aes(x=z, y=lambda)) +
    xlab("z") + ylab(expression(lambda)) +
    theme_bw(18) + theme(panel.grid=element_blank())
  ggsave(h1,file="~/Dropbox/Thesis/img/lambda_z1.png",
         height=4.5, width=4.5, units="in", dpi=300)

}



fret_example <- function(){
  #mu <-rep(c(0, 1.5, 0, 1.5, 0, 1.5, 0), c(50, 10, 50, 10, 80, 10, 40))
  mu <-rep(c(0, 1.5, 0, -1.5, 0), c(50, 25, 125, 10, 40))
  q <- rle(mu)
  ll <- cumsum(q$lengths)
  ll <- ll[-length(ll)]
  rect.data <- data.frame(x=ll,
                          type=rep(c("xmin", "xmax"), length(ll)/2),
                          w = rep(c(1:(length(ll)/2)), each=2))
  rect.data <- spread(rect.data, "type", "x")
  rect.data$xmax <- rect.data$xmax + 1
  rect.data$xmin <- rect.data$xmin
  rect.data$ymin <- -Inf
  rect.data$ymax=Inf
  set.seed(3*10^7)
  y <- rnorm(n=length(mu), mean=mu)
  ysm <- ksmooth(x=1:length(mu), y=y, x.points = 1:length(mu), bandwidth = 10)$y
  df <- data.frame(pos=1:length(mu), y=y, ysm=ysm)
  df$mu <- mu
  plot1 <- ggplot(df) + geom_rect(data = rect.data,
                                  aes(xmax=xmax, xmin=xmin, ymax=Inf, ymin=-Inf),
                                  fill="blue", alpha=0.2) +
          geom_point(aes(x=pos, y=y), shape=1) + geom_hline(yintercept = 0, linetype=3) +
        labs(x="Position", y="Statistic", title="Step 1: Calculate Statistics") +
        theme_bw(18) + theme(panel.grid=element_blank())
  ggsave(plot1, file="~/Dropbox/Thesis/img/example_fret1.png", height=5, width=6, units="in", dpi=300)


  plot1_talk <- ggplot(df) + geom_rect(data = rect.data,
                                  aes(xmax=xmax, xmin=xmin, ymax=Inf, ymin=-Inf),
                                  fill="blue", alpha=0.3, color="blue", lwd=0.8) +
    geom_point(aes(x=pos, y=y), shape=20, color="darkgrey", size=1.2) + geom_hline(yintercept = 0, linetype=3) +
    geom_point(aes(x=pos, y=mu), shape=20, size=0.7) +
    labs(x="Position", y="Statistic") +
    theme_bw(18) + theme(panel.grid=element_blank(),
                         panel.border=element_blank(), axis.line.x =element_line(),
                         axis.line.y=element_line(), plot.title=element_blank())

  ggsave(plot1_talk, file="~/Dropbox/Thesis/final_talk/img/example_fret1.png",
         height=4, width=7, units="in", dpi=300)


  df$type <- factor(rep(1, nrow(df)))
  plot2 <- ggplot(df) + geom_rect(data = rect.data,
                                  aes(xmax=xmax, xmin=xmin, ymax=Inf, ymin=-Inf),
                                  fill="blue", alpha=0.2) +
    geom_point(aes(x=pos, y=y, shape=type)) +
    scale_shape_manual(labels=c(expression("T("~s["j"]~")")), values=1) +
    geom_line(aes(x=pos, y=ysm, color=type), lwd=1) +
    scale_color_manual( labels=c(expression(tilde(T)~"("~s["j"]~")")), values=1) +
    geom_hline(yintercept = 0, linetype=3) +
    #annotate(geom="rect", ymin=-3.2, ymax=-1.2, xmin=0, xmax=25 , alpha=0, color="black", lwd=1)+
    labs(x="Position", y=expression("Statistic"), title="Step 2: Smooth") +
    theme_bw(18) + theme(panel.grid=element_blank(), legend.position=c(0.1, 0.85),
                         legend.title=element_blank(), legend.key=element_blank(), axis.title.y=element_blank())



  ggsave(plot2, file="~/Dropbox/Thesis/img/example_fret2.png", height=5, width=6, units="in", dpi=300)


  plot2_talk <- ggplot(df) + geom_rect(data = rect.data,
                                       aes(xmax=xmax, xmin=xmin, ymax=Inf, ymin=-Inf),
                                       fill="blue", alpha=0.3, color="blue", lwd=0.8) +
    geom_point(aes(x=pos, y=y), shape=20, color="darkgrey", size=1.2) + geom_hline(yintercept = 0, linetype=3) +
    #geom_point(aes(x=pos, y=mu), shape=20, size=0.7) +
    geom_line(aes(x=pos, y=ysm), lwd=1) +
    labs(x="Position", y="Statistic") +
    theme_bw(18) + theme(panel.grid=element_blank(),
                         panel.border=element_blank(), axis.line.x =element_line(),
                         axis.line.y=element_line(), plot.title=element_blank())

  ggsave(plot2_talk, file="~/Dropbox/Thesis/final_talk/img/example_fret2.png",
         height=4, width=7, units="in", dpi=300)


  #thresh <- 0.75
  thresh <- 0.85
  q <- rle(abs(df$ysm) >= thresh)
  ll <- cumsum(q$lengths)
  n <- sum(q$values==TRUE)
  rect.data2 <- data.frame(xmin=ll[q$values==FALSE][1:n]-1, xmax=ll[q$values==TRUE][1:n] + 1)

  plot3 <- ggplot(df) + geom_rect(data = rect.data,
                                  aes(xmax=xmax, xmin=xmin, ymax=Inf, ymin=-Inf),
                                  fill="blue", alpha=0.2) +
    geom_rect(data=rect.data2, aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf), fill="darkorange2", alpha=0.2) +
    geom_point(aes(x=pos, y=y), shape=1) +
    geom_line(aes(x=pos, y=ysm), lwd=1) +
    geom_hline(yintercept = 0, linetype=3) +
    geom_hline(yintercept = c(-1, 1)*thresh, linetype=2) +
    labs(x="Position", title="Step 3: Threshold") +
    theme_bw(18) + theme(panel.grid=element_blank(), axis.title.y=element_blank())
  ggsave(plot3, file="~/Dropbox/Thesis/img/example_fret3.png",
         height=5, width=6, units="in", dpi=300)


  plot3_talk <- ggplot(df) + geom_rect(data = rect.data,
                                       aes(xmax=xmax, xmin=xmin, ymax=Inf, ymin=-Inf),
                                       fill="blue", alpha=0.3, color="blue", lwd=0.8) +
    geom_rect(data=rect.data2, aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf),
              fill="violetRed", alpha=0.3, color="violetRed", lwd=0.8) +
    geom_point(aes(x=pos, y=y), shape=20, color="darkgrey", size=1.2) + geom_hline(yintercept = 0, linetype=3) +
    #geom_point(aes(x=pos, y=mu), shape=20, size=0.7) +
    geom_line(aes(x=pos, y=ysm), lwd=1) +
    geom_hline(yintercept = 0, linetype=3) +
    geom_hline(yintercept = c(-1, 1)*thresh, linetype=2) +
    annotate("text", label="z", x=250, y=thresh + 0.25, size=7) +
    labs(x="Position", y="Statistic") +
    theme_bw(18) + theme(panel.grid=element_blank(),
                         panel.border=element_blank(), axis.line.x =element_line(),
                         axis.line.y=element_line(), plot.title=element_blank())

  ggsave(plot3_talk, file="~/Dropbox/Thesis/final_talk/img/example_fret3.png",
         height=4, width=7, units="in", dpi=300)


}

plot_merge <- function(){
  mu <- rep(c(0, 1.3, 0), c(80, 80, 80))
  set.seed(3*10^7)
  y <- rnorm(n=length(mu), mean=mu)
  ysm <- ksmooth(x=1:length(mu), y=y, x.points = 1:length(mu), bandwidth = 10)$y
  df <- data.frame(pos=1:length(mu), y=y, ysm=ysm)

  thresh <- 1

  q <- rle(abs(df$ysm) >= thresh)
  ll <- cumsum(q$lengths)
  n <- sum(q$values==TRUE)
  seg.dat <- data.frame(x=ll[q$values==FALSE][1:n],
                        xend=ll[q$values==TRUE][1:n] )

  plot <- ggplot(df) + geom_line(aes(x=pos, y=ysm)) +
    geom_point(aes(x=pos, y=y), col="grey", size=0.8) +
    geom_hline(yintercept = thresh, linetype=2) +
    geom_hline(yintercept = 0.3*thresh, linetype=3) +
    geom_segment(data=seg.dat, aes(x=x, xend=xend, y=thresh, yend=thresh),
                 lwd=1.4, col="violetRed") + #ylim(range(ysm)) +
    labs(x="Position", y="Statistic") +
         #y=expression(tilde(T)~"("~s[j]~")")) +
    theme_bw(18) + theme(panel.grid=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.ticks.x = element_blank()) +
    annotate("text", label="z", x=50, y=1.18, size=6) +
    annotate("text", label="z[0]", x=50, y=0.46, size=6, parse=TRUE)
  ggsave(plot, file="~/Dropbox/Thesis/img/merge2.png", height=4.5, width=15*0.6, units="in", dpi=300)

}


var_z_example <- function(){

  mu <- rep(c(0, .75, 0, 0, 1.5, 0), c(95, 10, 95, 95, 10, 95))
  set.seed(4*10^7)
  y <- rnorm(n=400, mean=mu, sd=rep(c(0.25, 1), each=200))
  ysm <- ksmooth(x=1:400, y=y, bandwidth = 10, x.points = 1:400)$y
  #plot(y, cex=0.5)

  df <- data.frame(pos=1:400, y=y, ysm=ysm)

  line.dat <- data.frame(x=c(0, 0, 200, 200),
                         xend=c(200, 200, 400, 400),
                         y=c(0.5, -0.5, 1, -1),
                         yend=c(0.5, -0.5, 1, -1))

  rect.dat <- data.frame(xmax=c(106, 306), xmin=c(94, 294), ymax=c(Inf, Inf), ymin=c(-Inf, -Inf))

 h <- ggplot(df)  +
   geom_rect(data=rect.dat, aes(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin), fill="blue", alpha=0.2) +
   geom_point(aes(x=pos, y=y), shape=1) + geom_line(aes(x=pos, y=ysm), lwd=1.2) +
   ylab("Statistic") +  xlab("Position") + #geom_vline(xintercept = 200, linetype=3) +
   geom_segment(data=line.dat, aes(x=x, xend=xend, y=y, yend=yend), linetype=2) +
   theme_bw(18) + theme(panel.grid=element_blank())

  ggsave(h, file="~/Dropbox/Thesis/img/var_z_example.png",
         height=4.5, width=15, units="in", dpi=300)

  h_talk <- ggplot(df)  +
    geom_rect(data=rect.dat, aes(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin), fill="blue", alpha=0.2) +
    geom_point(aes(x=pos, y=y), shape=20, col="darkgrey", size=0.8) +
    geom_line(aes(x=pos, y=ysm), lwd=0.8) +
    ylab("Statistic") +  xlab("Position") + #geom_vline(xintercept = 200, linetype=3) +
    ylim(c(-2.5, max(df$y))) +
    geom_segment(data=line.dat, aes(x=x, xend=xend, y=y, yend=yend), linetype=2) +
    theme_bw(18) + theme(panel.grid=element_blank(), panel.border = element_blank(),
                         axis.line.x=element_line(),
                         axis.line.y=element_line())

  ggsave(h_talk, file="~/Dropbox/Thesis/final_talk/img/var_z_example.png",
         height=4, width=7, units="in", dpi=300)

  h_talk2 <- h_talk +
              geom_segment(x=0, xend=200, y=-1.8, yend=-1.8, col="blue", lwd=1.5) +
              geom_segment(x=200, xend=400, y=-1.8, yend=-1.8, col="chartreuse3", lwd=1.5) +
              geom_segment(x=0, xend=400, y=-2.1, yend=-2.1, col="black", lwd=1.5) +
              annotate(geom="text", x=130, y=-1.5, label="E[1]", parse=TRUE, col="blue", size=5) +
              annotate(geom="text", x=270, y=-1.5, label="E[2]", parse=TRUE, col="chartreuse3", size=5) +
              annotate(geom="text", x=200, y=-2.4, label="D", parse=TRUE, size=5)
  ggsave(h_talk2, file="~/Dropbox/Thesis/final_talk/img/var_z_example2.png",
         height=4, width=7, units="in", dpi=300)

}




region_cartoon <- function(){
  mu <- rep(c(0, 1.5, 0), c(7, 5, 8))
  df <- data.frame(pos=1:length(mu), mu=mu)
  p <- ggplot(df) +  geom_hline(yintercept = 0, lty=2) +
    geom_point(aes(x=pos, y=mu), shape=20, size=6, col="blue") +
    scale_y_continuous(limits=c(-2, 2.5)) +
    theme_bw() + theme(panel.grid=element_blank(), axis.text = element_blank(),
                       axis.title=element_blank(), axis.ticks=element_blank(),
                      panel.border=element_blank(), plot.title=element_blank())
  ggsave(p, file="~/Dropbox/Thesis/final_talk/img/betas.png",
         height=4, width=7, units="in", dpi=300)

  p1 <- p + geom_rect(xmin=7.5, xmax=12.5, ymin=-0.5, ymax=2, color="chartreuse3", alpha=0, lwd=1.5)
  ggsave(p1, file="~/Dropbox/Thesis/final_talk/img/betas1.png",
         height=4, width=7, units="in", dpi=300)


  p2 <- p + geom_rect(xmin=2.5, xmax=15.5, ymin=-0.5, ymax=2,
                      color="red", alpha=0, lwd=1.5)
  ggsave(p2, file="~/Dropbox/Thesis/final_talk/img/betas2.png",
         height=4, width=7, units="in", dpi=300)
  p3 <- p + geom_rect(xmin=9.5, xmax=11.5, ymin=-0.5, ymax=2,
                      color="red", alpha=0, lwd=1.5)
  ggsave(p3, file="~/Dropbox/Thesis/final_talk/img/betas3.png",
         height=4, width=7, units="in", dpi=300)
  p4 <- p + geom_rect(xmin=10.5, xmax=16.5, ymin=-0.5, ymax=2,
                      color="red", alpha=0, lwd=1.5)
  ggsave(p4, file="~/Dropbox/Thesis/final_talk/img/betas4.png",
         height=4, width=7, units="in", dpi=300)


  mu <- rep(c(0, 1, -0.5, 0), c(7, 5, 3, 8))
  df <- data.frame(pos=1:length(mu), mu=mu)
  p5 <- ggplot(df) +  geom_hline(yintercept = 0, lty=2) +
    geom_point(aes(x=pos, y=mu), shape=20, size=6, col="blue") +
    scale_y_continuous(limits=c(-2, 2.5)) +
    theme_bw() + theme(panel.grid=element_blank(), axis.text = element_blank(),
                       axis.title=element_blank(), axis.ticks=element_blank(),
                       panel.border=element_blank(), plot.title=element_blank())

  p5 <- p5 + geom_rect(xmin=7.5, xmax=15.5, ymin=-1, ymax=1.5, color="red",
                       alpha=0, lwd=1.5)
  ggsave(p5, file="~/Dropbox/Thesis/final_talk/img/betas5.png",
         height=4, width=7, units="in", dpi=300)

}

false_discovery_cartoon <- function(){

  mu <- rep(c(0, 1.5, 0, 1.5, 0), c(7, 5, 7, 8, 4))
  df <- data.frame(pos=1:length(mu), mu=mu)

  rect.data <- data.frame(xmin=c(7.5,19.5), xmax=c(12.5, 27.5 ),
                          ymin=rep(-1, 2), ymax=rep(2, 2))
  rect.data2 <- data.frame(xmin=c(1.5, 9.5, 21.5), xmax=c(3.5, 14.5, 24.5),
                           ymin=rep(-0.9, 3), ymax=rep(1.9, 3))

  p <-  ggplot(df) +  geom_rect(data=rect.data, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                                fill="blue", alpha=0.3, color="blue", lwd=1.5) +
    geom_hline(yintercept = 0, lty=2) +
    geom_point(aes(x=pos, y=mu), shape=20, size=6, col="blue") +
    scale_y_continuous(limits=c(-1, 2.5)) +
    theme_bw() + theme(panel.grid=element_blank(), axis.text = element_blank(),
                       axis.title=element_blank(), axis.ticks=element_blank(),
                       panel.border=element_blank(), plot.title=element_blank())

  ggsave(p, file="~/Dropbox/Thesis/final_talk/img/rfdr0.png",
         height=4, width=7, units="in", dpi=300)



  p1 <- ggplot(df) +  geom_rect(data=rect.data, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                               fill="blue", alpha=0.3, color="blue", lwd=1.5) +
    geom_rect(data=rect.data2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill="violetRed", alpha=0.3, color="violetRed", lwd=1.5) +
    geom_hline(yintercept = 0, lty=2) +
    geom_point(aes(x=pos, y=mu), shape=20, size=6, col="blue") +
    scale_y_continuous(limits=c(-1, 2.5)) +
    theme_bw() + theme(panel.grid=element_blank(), axis.text = element_blank(),
                       axis.title=element_blank(), axis.ticks=element_blank(),
                       panel.border=element_blank(), plot.title=element_blank())

  ggsave(p1, file="~/Dropbox/Thesis/final_talk/img/rfdr1.png",
         height=4, width=7, units="in", dpi=300)


  rect.data2 <- data.frame(xmin=c(3.5), xmax=c(28.5),
                           ymin=rep(-0.9, 1), ymax=rep(1.9, 1))

  p2 <- ggplot(df) +  geom_rect(data=rect.data, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                               fill="blue", alpha=0.3, col="blue", lwd=1.5) +
    geom_rect(data=rect.data2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill="violetRed", alpha=0.3, col="violetRed", lwd=1.5) +
    geom_hline(yintercept = 0, lty=2) +
    geom_point(aes(x=pos, y=mu), shape=20, size=6, col="blue") +
    scale_y_continuous(limits=c(-1, 2.5)) +
    theme_bw() + theme(panel.grid=element_blank(), axis.text = element_blank(),
                       axis.title=element_blank(), axis.ticks=element_blank(),
                       panel.border=element_blank(), plot.title=element_blank())

  ggsave(p2, file="~/Dropbox/Thesis/final_talk/img/rfdr2.png",
         height=4, width=7, units="in", dpi=300)


  mu <- rep(c(0, 0, 0, 1.5, 0), c(7, 5, 7, 10, 2))
  df <- data.frame(pos=1:length(mu), mu=mu)

  rect.data <- data.frame(xmin=c(19.5), xmax=c( 29.5 ),
                          ymin=rep(-1, 1), ymax=rep(2, 1))
  rect.data2 <- data.frame(xmin=c(1.5, 20.5, 22.5, 24.5, 26.5, 28.5),
                           xmax=c(17.5, 21.5, 23.5, 25.5, 27.5, 29.5),
                           ymin=rep(-0.9, 6), ymax=rep(1.9, 6))

  p3 <- ggplot(df) +  geom_rect(data=rect.data, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                                fill="blue", alpha=0.3, col="blue", lwd=1.5) +
    geom_rect(data=rect.data2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill="violetRed", alpha=0.3, col="violetRed", lwd=1.5) +
    geom_hline(yintercept = 0, lty=2) +
    geom_point(aes(x=pos, y=mu), shape=20, size=6, col="blue") +
    scale_y_continuous(limits=c(-1, 2.5)) +
    theme_bw() + theme(panel.grid=element_blank(), axis.text = element_blank(),
                       axis.title=element_blank(), axis.ticks=element_blank(),
                       panel.border=element_blank(), plot.title=element_blank())

  ggsave(p3, file="~/Dropbox/Thesis/final_talk/img/rfdr3.png",
         height=4, width=7, units="in", dpi=300)

}


plot_nregion <- function(){

  #FRET
  df = getobj("region_counts.RData")
  #names(df)= c("fdr", "total",  "Regions overlapping peaks")
  df$`All regions` = df$total - df$peak
  df$`Regions overlapping peaks` <- df$peak - df$contained
  df$`Regions contained in peaks` <- df$contained
  ymax = max(df$total)
  df <- df[, c("fdr", "All regions", "Regions overlapping peaks", "Regions contained in peaks")]

  dflong = gather(df, "class", "count", -fdr )
  #dflong$class=factor(dflong$class, levels=c("All regions", "Regions overlapping peaks", "Regions contained in peaks")   )
  dflong$class=factor(dflong$class, levels=c("Regions contained in peaks", "Regions overlapping peaks", "All regions")   )
  o <- c(241:360, 121:240, 1:120)
  dflong <- dflong[o, ]
  gridlines = seq(10000, ymax, by=10000)
  h_fret = ggplot(dflong) + geom_area(aes(x=fdr, y = count, fill=class), alpha=0.5) +
    ylab("Number of Regions") + xlab("rFDR Threshold") +
    scale_x_continuous(limits=c(0.03, 0.2))+
    scale_fill_manual(values=c("chartreuse3", "blue", "darkorange2"))+
    geom_hline(yintercept = gridlines, lty=3) +ggtitle("FRET") +
    theme_bw(14) + theme(panel.grid=element_blank(),
                         legend.text = element_text(size=11),
                         legend.background = element_blank(),
                         legend.position=c(0.35, 0.9), legend.title=element_blank())

  ggsave(h_fret, file="~/Dropbox/Thesis/final_talk/img/fret_results2.png",
         height=4.5, width=4.5, units="in", dpi=300)

  #Wellington Bootstrap
  df = getobj("region_counts_well.RData")
  names(df)= c("score", "total",  "Footprints overlapping peaks")
  df=df[df$score >=50,]
  ymax = max(df$total)
  df$`All footprints` = df$total - df$`Footprints overlapping peaks`

  dflong = gather(df, "class", "count", -score ,-total)
  dflong$class=factor(dflong$class, levels=rev(levels(dflong$class)))
  gridlines = seq(10000, ymax, by=10000)
  h_well = ggplot(dflong) + geom_area(aes(x=score, y = count, fill=class), alpha=0.5) +
    ylab("Number of Footprints") + xlab("Score") +
    scale_x_reverse(breaks=seq(10, 100, by=20), limits=c(100, 50))+
    scale_fill_manual(values=c("darkorange1", "blue"))+
    geom_hline(yintercept = gridlines, lty=3) +
    theme_bw(12) + theme(panel.grid=element_blank(),
                         legend.position=c(0.3, 0.9), legend.title=element_blank())

  #Huber
  df = getobj("region_counts_huber.RData")
  names(df)= c("fdr", "total")
  ymax = max(df$total)
  gridlines = seq(10000, ymax, by=10000)
  h_huber = ggplot(df) + geom_area(aes(x=fdr, y = total), fill="blue", alpha=0.5) +
    ylab("Number of Peaks") + xlab("FDR Threshold") +
    geom_hline(yintercept = gridlines, lty=3) +
    scale_x_continuous(breaks=seq(0, 0.2, by=0.02), limits=c(0, 0.2))+
    theme_bw(12) + theme(panel.grid=element_blank())

  h_huber_talk = ggplot(df) + geom_area(aes(x=fdr, y = total), fill="blue", alpha=0.5) +
    ylab("Number of Peaks") + xlab("FDR Threshold") +
    geom_hline(yintercept = gridlines, lty=3) +
    scale_x_continuous(limits=c(0, 0.2))+
    theme_bw(14) + theme(panel.grid=element_blank()) + ggtitle("Huber")
  ggsave(h_huber_talk, file="~/Dropbox/Thesis/final_talk/img/huber_results.png",
         height=4.5, width=4.5, units="in", dpi=300)


  #DESeq2
  df = getobj("region_counts_deseq2.RData")
  names(df)= c("fdr", "total")
  ymax = max(df$total)
  gridlines = seq(10000, ymax, by=10000)
  h_deseq2 = ggplot(df) + geom_area(aes(x=fdr, y = total), fill="blue", alpha=0.5) +
    ylab("Number of Peaks") + xlab("FDR Threshold") +
    scale_x_continuous(breaks=seq(0, 0.2, by=0.02), limits=c(0, 0.2))+
    geom_hline(yintercept = gridlines, lty=3) +
    theme_bw(12) + theme(panel.grid=element_blank())
  h_deseq2_talk = ggplot(df) + geom_area(aes(x=fdr, y = total), fill="blue", alpha=0.5) +
    ylab("Number of Peaks") + xlab("FDR Threshold") +
    scale_x_continuous(limits=c(0, 0.2))+
    geom_hline(yintercept = gridlines, lty=3) +
    theme_bw(14) + theme(panel.grid=element_blank()) + ggtitle("DESeq2")
  ggsave(h_deseq2_talk, file="~/Dropbox/Thesis/final_talk/img/deseq2_results.png",
         height=4.5, width=4.5, units="in", dpi=300)

  #Waveqtl
  df = getobj("region_counts_wave.RData")
  names(df)= c("pval", "total")
  ymax = max(df$total)
  gridlines = seq(10000, ymax, by=10000)
  h_wave = ggplot(df) + geom_area(aes(x=pval, y = total), fill="blue", alpha=0.5) +
    ylab("Number of Peaks") + xlab("p-value Threshold") +
    geom_hline(yintercept = gridlines, lty=3) +
    theme_bw(12) + theme(panel.grid=element_blank())


  #Titles and scales
  h_fret = h_fret + ggtitle("FRET")
  h_well = h_well + ggtitle("Wellington-Bootstrap")
  h_huber = h_huber + ggtitle("Huber Fixed-Window Test")
  h_deseq2 = h_deseq2 + ggtitle("DESeq2")
  h_wave = h_wave + ggtitle("WaveQTL")

  ggsave(h_fret, file="~/Dropbox/Thesis/img/fret_results2.png", height=4.5, width=4.5, units="in", dpi=300)

  ggsave(h_well, file="~/Dropbox/Thesis/img/well_results.png", height=4.5, width=4.5, units="in", dpi=300)
  ggsave(h_huber, file="~/Dropbox/Thesis/img/huber_results.png", height=4.5, width=4.5, units="in", dpi=300)
  ggsave(h_deseq2, file="~/Dropbox/Thesis/img/deseq2_results.png", height=4.5, width=4.5, units="in", dpi=300)
  ggsave(h_wave, file="~/Dropbox/Thesis/img/wave_results.png", height=4.5, width=4.5, units="in", dpi=300)
}


sim_results_for_talk <- function(){
  tot.rates <- getobj("sim_results/bin1_total2.RData")

  qqs <- list()
  for(i in 1:20) qqs[[i]] <- cfdrSims:::thin_qqplot(pvals=tot.rates$null.p[[i+9]]) +
    theme(text=element_text(size=18)) + ggtitle(tot.rates$names[i+9])

  h <- grid.arrange(qqs[[1]], qqs[[2]], qqs[[3]], qqs[[4]], qqs[[5]],
                    qqs[[6]], qqs[[7]], qqs[[8]], qqs[[9]], qqs[[10]],
                    qqs[[11]], qqs[[12]], qqs[[13]], qqs[[14]], qqs[[15]],
                    qqs[[16]], qqs[[17]], qqs[[18]], qqs[[19]], qqs[[20]], ncol=5)
  ggsave(h, file="bin2_qqs.png", heigh=10, width=10, units="in", dpi=300)

  nn <- "w64bDESeq2"
  ix <- which(tot.rates$names==nn)
  h1 <- cfdrSims:::thin_qqplot(pvals=tot.rates$null.p[[ix]]) +
    theme(text=element_text(size=12))
  ggsave(h1, file="~/Dropbox/Thesis/final_talk/img/deseq2_64b_qq.png", height=3, width=3, units="in")

  nn <- "w64eDESeq2"
  ix <- which(tot.rates$names==nn)
  h2<- cfdrSims:::thin_qqplot(pvals=tot.rates$null.p[[ix]]) +
    theme(text=element_text(size=12))
  ggsave(h2, file="~/Dropbox/Thesis/final_talk/img/deseq2_64e_qq.png", height=3, width=3, units="in")

  p1 <- make_sim_plot(tot.rates, names=c("fretHuber2", "fretHuber1", "fretHuber6"),
                      col=c("black", "violetRed", "chartreuse3"),
                      lty=c(1, 1, 1), shapes=c(15, 15, 15))
  ggsave(p1$tprplot, file="~/Dropbox/Thesis/final_talk/img/fret_compare2_tpp.png",height=4, width=4, units="in", dpi=300)
  ggsave(p1$fdpplot, file="~/Dropbox/Thesis/final_talk/img/fret_compare2_fdp.png",height=4, width=4, units="in", dpi=300)

  p2 <- make_sim_plot(tot.rates, names=c("fretHuber2", "fretPois2", "fretT2"),
                      lty=c(1, 2, 3), col=rep("black", 3), shapes=c(15, 16, 17))
  ggsave(p2$tprplot, file="~/Dropbox/Thesis/final_talk/img/fret_compare_tpp.png",height=4, width=4, units="in", dpi=300)
  ggsave(p2$fdpplot, file="~/Dropbox/Thesis/final_talk/img/fret_compare_fdp.png",height=4, width=4, units="in", dpi=300)

  l <- make_sim_legend_talk()
  ggsave(l$plot, file="~/Dropbox/Thesis/final_talk/img/legend.png",
         height=4, width=3, units="in", dpi=300)

  p64b <- make_sim_plot(tot.rates, names=l$info$informed,
                        cols=l$info$cols,
                        ltys=l$info$lty, shapes=l$info$shape)
  ggsave(p64b$fdpplot, file="~/Dropbox/Thesis/final_talk/img/w64b_fdp.png", height=4, width=4, units="in", dpi=300)
  ggsave(p64b$tprplot, file="~/Dropbox/Thesis/final_talk/img/w64b_tpp.png", height=4, width=4, units="in", dpi=300)


  p64e <- make_sim_plot(tot.rates, names=l$info$naive,
                        cols=l$info$cols,
                        ltys=l$info$lty, shapes=l$info$shape)
  ggsave(p64e$fdpplot, file="~/Dropbox/Thesis/final_talk/img/w64e_fdp.png", height=4, width=4, units="in", dpi=300)
  ggsave(p64e$tprplot, file="~/Dropbox/Thesis/final_talk/img/w64e_tpp.png", height=4, width=4, units="in", dpi=300)


  p64o <- make_sim_plot(tot.rates, names=c("fretHuber2", "w64oDESeq2", "w64oWave", "w64oHuber", "w64oPois"),
                        cols=c("black", "blue", "violetRed", "chartreuse3", "darkorange2"),
                        ltys=rep(1, 5), shapes=15:19)

  p8b <- make_sim_plot(tot.rates, names=c("fretHuber2", "w8bDESeq2", "w8bWave", "w8bHuber", "w8bPois"),
                        cols=c("black", "blue", "violetRed", "chartreuse3", "darkorange2"),
                        ltys=rep(1, 5), shapes=15:19)


  pdeseq2 <- make_sim_plot(tot.rates, names=c("w64bDESeq2", "w64eDESeq2", "w64oDESeq2", "w8bDESeq2"),
                                               cols=rep("blue", 4),
                                               ltys=1:4, shapes=15:18)

  pwave <- make_sim_plot(tot.rates, names=c("w64bWave", "w64eWave", "w64oWave", "w8bWave"),
                           cols=rep("violetRed", 4),
                           ltys=1:4, shapes=15:18)

  phuber <- make_sim_plot(tot.rates, names=c("w64bHuber", "w64eHuber", "w64oHuber", "w8bHuber"),
                         cols=rep("chartreuse3", 4),
                         ltys=1:4, shapes=15:18)

}



make_sim_legend_talk <- function(){
  points <- data.frame(x=rep(1, 5), y = rev(seq(1, 3, length.out=5)),
                       x2 = 1:5, y2=rep(1, 5))

  points$left = points$x2-0.45
  points$right = points$x2 + 0.45

  points$labs <- c("FRET", "DESeq2", "WaveQTL", "Huber", "QP")
  points$lty <- c(1, 2, 4, 2, 3)
  points$shape <- 15:19
  points$informed <- c("fretHuber2", paste0("w64b", c("DESeq2", "Wave", "Huber", "Pois")))
  points$naive <- c("fretHuber2", paste0("w64e", c("DESeq2", "Wave", "Huber", "Pois")))

  points$cols <- c("black", "blue", "violetRed", "chartreuse3", "darkorange2")

  p <-ggplot(points) +
    geom_segment(aes(x=left, xend=right, y=y2, yend=y2), lty=points$lty,
                 colour=points$cols, lwd=1.5)+
    geom_point(aes(x=x2, y=y2), col="white", shape=20, size=13)+
    geom_point(aes(x=x2, y=y2), col=points$cols, shape=points$shape, size=2, stroke=1.3)+
    annotate(geom="text", y=rep(1.001, 5), x=points$x2,
             label=points$labs, size=5) +
    xlim(0.55, 5.45) +
    theme_bw() + theme(panel.grid=element_blank(), axis.title=element_blank(),
                       axis.text=element_blank(),
                       panel.border=element_blank(), axis.ticks=element_blank())

  #ggsave(p, file="~/Dropbox/Thesis/final_talk/img/legend_h.png", height=2, width=7, units="in", dpi=300)
  return(list("plot"=p, "info"=points))
}

make_sim_legend_talk2 <- function(){
  points <- data.frame(y=rep(1, 3), x = 1:3)

  points$left = points$x-0.45
  points$right = points$x + 0.45

  #points$labs <- c("FRET-Huber", "FRET-t-test", "FRET-QP")
  points$labs <- c("FRET-1", "FRET-2", "FRET-6")
  #points$lty <- c(1, 2, 3)
  points$lty <- c(1, 1, 1)
  #points$shape <- 15:17
  points$shape <- 15
  points$names <- c("fretHuber1", "fretHuber2", "fretHuber6")

  points$cols <- c("violetRed", "black", "chartreuse3")

  p <-ggplot(points) +
    geom_segment(aes(x=left, xend=right, y=y, yend=y), lty=points$lty,
                 colour=points$cols, lwd=1.5)+
    geom_point(aes(x=x, y=y), col="white", shape=20, size=13)+
    geom_point(aes(x=x, y=y), col=points$cols, shape=points$shape, size=2, stroke=1.3)+
    annotate(geom="text", y=rep(1.001, 3), x=points$x,
             label=points$labs, size=5) +
    xlim(0.55, 3.45) +
    theme_bw() + theme(panel.grid=element_blank(), axis.title=element_blank(),
                       axis.text=element_blank(),
                       panel.border=element_blank(), axis.ticks=element_blank())

  #ggsave(p, file="~/Dropbox/Thesis/final_talk/img/legend_h_Nfret.png", height=2, width=7, units="in", dpi=300)
  return(list("plot"=p, "info"=points))
}

loss_plot <- function(){
  x <- seq(-2, 2, length.out=200)
  cols <- c("darkorange", "chartreuse3", "blue")
  h <- function(x){
    if(abs(x) < 1) return(0.5*x^2)
    return(abs(x)-0.5)
  }
  df <- data.frame(x=x, yse= 0.5*x^2, yme=abs(x), yhu = sapply(x, FUN=h) )
  dflong <- gather(df, "Loss", "y", -x)
  p <- ggplot(dflong) + geom_line(aes(x=x, y=y, group=Loss, color=Loss, linetype=Loss), lwd=2) +
    scale_color_manual(labels=c("Squared Error", "Median", "Huber, c=1"), values=cols) +
    scale_linetype_manual(labels=c("Squared Error", "Median", "Huber, c=1"), values=c(1, 4, 2))+
    xlab(expression(u)) + ylab(expression(paste(rho, "(", u, ")"))) +
    theme_bw(14) + theme(panel.grid=element_blank(),
                         legend.title=element_blank(),
                         legend.key=element_blank())
  ggsave(p, file="~/Dropbox/Thesis/final_talk/img/loss.png", height=4.5, width=6, units="in", dpi=300)
}


st_ts_roc_curves_normal <- function(){

  #Normal
  R20 <- getobj("~/Desktop/smooth_test/ar_sd2_rho_0_bw20_ts_ttest.RData")
  R5 <- getobj("~/Desktop/smooth_test/ar_sd2_rho_0_bw5_ts_ttest.RData")
  R50 <- getobj("~/Desktop/smooth_test/ar_sd2_rho_0_bw50_ts_ttest.RData")

  #which.ci <- seq(0, 200, length.out=6)

  #plotCI(x=my.interp$fpr[whichCI], y=my.interp$tpr[whichCI],
  #       uiw=my.interp$s.e[whichCI], err=err, col=c, add=TRUE, pch=0)
  #lines(x=my.interp$fpr, y=my.interp$tpr, col=c, lwd=lwd, lty=l)


  pw_df_normal<- data.frame(fpr = rep(R20$st_rates[[1]]$fpr, 6),
                          bw = rep(c(20, 5, 50), each=400),
                         type = rep(  rep(c("st", "ts"), each=200) , 3),
                         rates = c( R20$st_rates[[1]]$tpr, R20$ts_rates[[1]]$tpr,
                                    R5$st_rates[[1]]$tpr, R5$ts_rates[[1]]$tpr,
                                    R50$st_rates[[1]]$tpr, R50$ts_rates[[1]]$tpr),
                         se = c(R20$st_rates[[1]]$s.e, R20$ts_rates[[1]]$s.e,
                                R5$st_rates[[1]]$s.e, R5$ts_rates[[1]]$s.e,
                                R50$st_rates[[1]]$s.e, R50$ts_rates[[1]]$s.e))
  #which.ci = ceiling(seq(1, 200, length.out=7)[2:6]
  pw_df_normal$bw <- as.factor(pw_df_normal$bw)
  pw_normal <- ggplot(pw_df_normal) + geom_line(aes(x=fpr, y=rates, col=bw, lty=type), lwd=1) +
    #geom_errorbar(data=pw_df_normal[which.ci,], aes(x=fpr, ymax=rates + se, ymin=rates-se, col=bw)) +
    xlab("False Positive Rate") + ylab("Avg. True Detection Proportion") +
    scale_color_manual(values=c("blue", "black", "darkorange2")) +
    scale_linetype_manual(values=c(2, 1)) + theme_bw(16) + theme(panel.grid=element_blank(),
                                                                 legend.position="none")

  ggsave(pw_normal, file="~/Dropbox/Thesis/img/st_normal_pw.png", height=4, width=4, units="in", dpi=300)

  n20_st <- length(R20$st_rates[[2]]$fct)
  n20_ts <- length(R20$ts_rates[[2]]$fct)
  n20 <- n20_ts + n20_st
  n50_st <- length(R50$st_rates[[2]]$fct)
  n50_ts <- length(R50$ts_rates[[2]]$fct)
  n50 <- n50_ts + n50_st
  n5_st <- length(R5$st_rates[[2]]$fct)
  n5_ts <- length(R5$ts_rates[[2]]$fct)
  n5 <- n5_ts + n5_st
  rw_df_normal<- data.frame(fct = c( R20$st_rates[[2]]$fct, R20$ts_rates[[2]]$fct,
                                     R5$st_rates[[2]]$fct, R5$ts_rates[[2]]$fct,
                                     R50$st_rates[[2]]$fct, R50$ts_rates[[2]]$fct),
                            bw = rep(c(20, 5, 50), c(n20, n5, n50)),
                            type = rep(rep(c("st", "ts"), 3), c(n20_st, n20_ts, n5_st, n5_ts, n50_st, n50_ts)),
                            rates = c( R20$st_rates[[2]]$tpr, R20$ts_rates[[2]]$tpr,
                                       R5$st_rates[[2]]$tpr, R5$ts_rates[[2]]$tpr,
                                       R50$st_rates[[2]]$tpr, R50$ts_rates[[2]]$tpr))
  rw_df_normal$bw <- as.factor(rw_df_normal$bw)
  rw_normal <- ggplot(rw_df_normal) + geom_line(aes(x=fct, y=rates, col=bw, lty=type), lwd=1) +
    #geom_errorbar(data=pw_df_normal[which.ci,], aes(x=fpr, ymax=rates + se, ymin=rates-se, col=bw)) +
    xlab("Number of False Discoveries") + ylab("Avg. Prop. of Regions Detected") +
    scale_color_manual(values=c("blue", "black", "darkorange2")) +
    scale_linetype_manual(values=c(2, 1)) + theme_bw(16) + theme(panel.grid=element_blank(),
                                                                 legend.position="none")
  ggsave(rw_normal, file="~/Dropbox/Thesis/img/st_normal_rw.png", height=4, width=4, units="in", dpi=300)

}



st_ts_roc_curves_poisson <- function(){

  #Poisson
  R20 <- getobj("~/Desktop/smooth_test/bin1_bw20_ts_t.RData")
  R5 <- getobj("~/Desktop/smooth_test/bin1_bw5_ts_t.RData")
  R50 <- getobj("~/Desktop/smooth_test/bin1_bw50_ts_t.RData")

  pw_df_normal<- data.frame(fpr = rep(R20$st_rates[[1]]$fpr, 6),
                            bw = rep(c(20, 5, 50), each=400),
                            type = rep(  rep(c("st", "ts"), each=200) , 3),
                            rates = c( R20$st_rates[[1]]$tpr, R20$ts_rates[[1]]$tpr,
                                       R5$st_rates[[1]]$tpr, R5$ts_rates[[1]]$tpr,
                                       R50$st_rates[[1]]$tpr, R50$ts_rates[[1]]$tpr))
  pw_df_normal$bw <- as.factor(pw_df_normal$bw)
  pw_normal <- ggplot(pw_df_normal) + geom_line(aes(x=fpr, y=rates, col=bw, lty=type), lwd=1) +
    #geom_errorbar(data=pw_df_normal[which.ci,], aes(x=fpr, ymax=rates + se, ymin=rates-se, col=bw)) +
    xlab("False Positive Rate") + ylab("Avg. True Detection Proportion") +
    scale_color_manual(values=c("blue", "black", "darkorange2")) +
    scale_linetype_manual(values=c(2, 1)) + theme_bw(16) + theme(panel.grid=element_blank(),
                                                                 legend.position="none")

  ggsave(pw_normal, file="~/Dropbox/Thesis/img/st_poisson_pw.png", height=4, width=4, units="in", dpi=300)

  n20_st <- length(R20$st_rates[[2]]$fct)
  n20_ts <- length(R20$ts_rates[[2]]$fct)
  n20 <- n20_ts + n20_st
  n50_st <- length(R50$st_rates[[2]]$fct)
  n50_ts <- length(R50$ts_rates[[2]]$fct)
  n50 <- n50_ts + n50_st
  n5_st <- length(R5$st_rates[[2]]$fct)
  n5_ts <- length(R5$ts_rates[[2]]$fct)
  n5 <- n5_ts + n5_st
  rw_df_normal<- data.frame(fct = c( R20$st_rates[[2]]$fct, R20$ts_rates[[2]]$fct,
                                     R5$st_rates[[2]]$fct, R5$ts_rates[[2]]$fct,
                                     R50$st_rates[[2]]$fct, R50$ts_rates[[2]]$fct),
                            bw = rep(c(20, 5, 50), c(n20, n5, n50)),
                            type = rep(rep(c("st", "ts"), 3), c(n20_st, n20_ts, n5_st, n5_ts, n50_st, n50_ts)),
                            rates = c( R20$st_rates[[2]]$tpr, R20$ts_rates[[2]]$tpr,
                                       R5$st_rates[[2]]$tpr, R5$ts_rates[[2]]$tpr,
                                       R50$st_rates[[2]]$tpr, R50$ts_rates[[2]]$tpr))
  rw_df_normal$bw <- as.factor(rw_df_normal$bw)
  rw_normal <- ggplot(rw_df_normal) + geom_line(aes(x=fct, y=rates, col=bw, lty=type), lwd=1) +
    #geom_errorbar(data=pw_df_normal[which.ci,], aes(x=fpr, ymax=rates + se, ymin=rates-se, col=bw)) +
    xlab("Number of False Discoveries") + ylab("Avg. Prop. of Regions Detected") +
    scale_color_manual(values=c("blue", "black", "darkorange2")) +
    scale_linetype_manual(values=c(2, 1)) + theme_bw(16) + theme(panel.grid=element_blank(),
                                                                 legend.position="none")
  ggsave(rw_normal, file="~/Dropbox/Thesis/img/st_poisson_rw.png", height=4, width=4, units="in", dpi=300)

}




make_st_legend <- function(){
  points <- cbind(rep(c(1, 2.2, 3.4), each=2), rep(c(1.03, 1.08), 3))
  points <- data.frame(points)
  names(points)=c("x", "y")
  points$left = points$x-0.45
  points$right = points$x + 0.65

  points$lty = 1
  points$lty[points$y==1.03] <- 2


  points$cols <- rep(c("blue", "black", "darkorange2"), each=2)

  p <-ggplot(points) +
    geom_segment(aes(x=left, xend=right, y=y, yend=y), lty=points$lty,
                 colour=points$cols, lwd=1)+
    annotate(geom="text", x=c(0, 1, 2.2, 3.4), y=rep(1.12, 4),
             label=c("Bandwidth:", "5", "20", "50"), size=4)+
    annotate(geom="text", x=c(4.4, 4.4), y=c(1.03, 1.08),
             label=c("ST", "TS"), size=4)+
    xlim(-0.5, 4.5) +
    theme_bw() + theme(panel.grid=element_blank(), axis.title=element_blank(),
                       axis.text=element_blank(),
                       panel.border=element_blank(), axis.ticks=element_blank())
  return(list(p, points))
}





