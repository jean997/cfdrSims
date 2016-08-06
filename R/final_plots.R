example_profile <- function(){
  yy <- cfdrSims:::gen_profile(5)
  yy <- c(yy,  cfdrSims:::gen_profile(3))
  yy <- c(yy,  cfdrSims:::gen_profile(7))
  yy <- c(yy,  cfdrSims:::gen_profile(c(3.2, 4.8)))
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

}


plot_lambda <- function(){
  R <- getobj("cont1_100_fret.RData")
  #signal <- Intervals(R$signal$signal)
  #perm.stats <- R$stats[2, , -1]
  #z1_true <- choose_z_even(perm.stats, nlam=50, bw=20, pos=1:12000,
                          #z0=R$cl[[2]][[1]]$z0, seg.ends=c(12000), except=signal)

  dfz1 <- data.frame(R$cl[[2]][[1]]$z)
  names(dfz1) = c("lambda", "z")
  dfz1$lambda <- 12000 * 10^dfz1$lambda
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


  L <- 0.6e-1
  z0 <- approx(x=dfz1$lambda, y=dfz1$z, xout=L)$y

  line.dat1 = data.frame(x = c(-Inf, z0, z0),  y = c(L, L, 0))

  h1_talk <- ggplot(dfz1) + geom_line(aes(x=z, y=lambda)) +
    geom_line(data=line.dat1, aes(x=x, y=y), lty=2) +
    geom_point(aes(x=z0, y=L), size=4) +
    xlab(expression("z")) + ylab(expression(lambda~"(z)")) +
    scale_y_log10(breaks=L, labels=expression(Lambda))+
    scale_x_continuous(breaks=z0, labels=expression(tilde("z"))) +
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






  labs <- c("Interval 1", "Interval 2")

  L <- 0.3e-1
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
    ylab(expression(lambda~"(z)")) +
    theme_bw(18) + theme(panel.grid=element_blank(),
                         legend.title=element_blank(),
                         panel.border = element_blank(),
                         axis.line.x=element_line(), #axis.text=element_blank(),
                         axis.line.y=element_line(),
                         legend.position="none")

  ggsave(h2_talk, file="~/Dropbox/Thesis/final_talk/img/lambda_2.png",
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
                                  fill="blue", alpha=0.2) +
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
                                       fill="blue", alpha=0.2) +
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
                                       fill="blue", alpha=0.2) +
    geom_rect(data=rect.data2, aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf), fill="darkorange2", alpha=0.2) +
    geom_point(aes(x=pos, y=y), shape=20, color="darkgrey", size=1.2) + geom_hline(yintercept = 0, linetype=3) +
    #geom_point(aes(x=pos, y=mu), shape=20, size=0.7) +
    geom_line(aes(x=pos, y=ysm), lwd=1) +
    geom_hline(yintercept = 0, linetype=3) +
    geom_hline(yintercept = c(-1, 1)*thresh, linetype=2) +
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
  rect.data2 <- data.frame(xmin=c(1.5, 9.5), xmax=c(3.5, 14.5),
                           ymin=rep(-0.9, 2), ymax=rep(1.9, 2))

  p <-  ggplot(df) +  geom_rect(data=rect.data, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                                fill="blue", alpha=0.3) +
    geom_hline(yintercept = 0, lty=2) +
    geom_point(aes(x=pos, y=mu), shape=20, size=6, col="blue") +
    scale_y_continuous(limits=c(-1, 2.5)) +
    theme_bw() + theme(panel.grid=element_blank(), axis.text = element_blank(),
                       axis.title=element_blank(), axis.ticks=element_blank(),
                       panel.border=element_blank(), plot.title=element_blank())

  ggsave(p, file="~/Dropbox/Thesis/final_talk/img/rfdr0.png",
         height=4, width=7, units="in", dpi=300)



  p1 <- ggplot(df) +  geom_rect(data=rect.data, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                               fill="blue", alpha=0.3) +
    geom_rect(data=rect.data2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill="orange", alpha=0.3) +
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
                               fill="blue", alpha=0.3) +
    geom_rect(data=rect.data2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill="orange", alpha=0.3) +
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
                                fill="blue", alpha=0.3) +
    geom_rect(data=rect.data2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill="orange", alpha=0.3) +
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
    scale_x_continuous(breaks=seq(0.04, 0.2, by=0.02), limits=c(0.03, 0.2))+
    scale_fill_manual(values=c("darkorange1", "blue", "chartreuse3"))+
    geom_hline(yintercept = gridlines, lty=3) +
    theme_bw(12) + theme(panel.grid=element_blank(),
                         legend.position=c(0.3, 0.9), legend.title=element_blank())

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
