
sim_plots_aoas <- function(){
  #tot.rates <- getobj("sim_results/bin1_total3.RData")

  qqs <- list()
  n <- length(tot.rates$nseg)*3
  for(i in (n + 1):length(tot.rates$names)) qqs[[i-n]] <- cfdrSims:::thin_qqplot(pvals=tot.rates$null.p[[i]]) +
    theme(text=element_text(size=18)) + ggtitle(tot.rates$names[i])

  h <- grid.arrange(qqs[[1]], qqs[[2]], qqs[[3]], qqs[[4]], qqs[[5]],
                    qqs[[6]], qqs[[7]], qqs[[8]], qqs[[9]], qqs[[10]], ncol=4)
  ggsave(h, file="bin1_qqs_64.png", heigh=10, width=10, units="in", dpi=300)

  #nn <- "w64bDESeq2"
  #ix <- which(tot.rates$names==nn)
  #h1 <- cfdrSims:::thin_qqplot(pvals=tot.rates$null.p[[ix]]) +
  #  theme(text=element_text(size=12))
  #ggsave(h1, file="~/Dropbox/Thesis/final_talk/img/deseq2_64b_qq.png", height=3, width=3, units="in")

  #nn <- "w64eDESeq2"
  #ix <- which(tot.rates$names==nn)
  #h2<- cfdrSims:::thin_qqplot(pvals=tot.rates$null.p[[ix]]) +
  #  theme(text=element_text(size=12))
  #ggsave(h2, file="~/Dropbox/Thesis/final_talk/img/deseq2_64e_qq.png", height=3, width=3, units="in")

  prefix <- "sim_results/bin4_nt2_upd"
  tot <- getobj(paste0(prefix, "_total.RData"))
  p1 <- make_sim_plot(tot, names=c("fretHuber2", "fretHuber1", "fretHuber20", "fretHuber60"),
                      col=c("black", "violetRed", "chartreuse3", "blue"),
                      lty=rep(1, 4), shapes=rep(15, 4))
  h <- grid.arrange(p1$fdpplot, p1$tprplot, ncol=2)
  ggsave(h, file=paste0(prefix, "_nseg.png"), height=4, width=8, units="in", dpi=300)

  #[1] "fretHuberauto50"  "fretHuberauto100" "fretHuberauto150"
  #[4] "fretHuberauto200" "fretHuberauto300" "fretHuberauto400"
  p1 <- make_sim_plot(tot, names=paste0("fretHuberauto", c(50, 100, 150, 200, 300, 400)),
                      col=c("black", "violetRed", "chartreuse3", "blue", "purple", "darkorange2"),
                      lty=c(1, 2, 1, 2, 1, 2), shapes=rep(15, 6))
  h <- grid.arrange(p1$fdpplot, p1$tprplot, ncol=2)
  ggsave(h, file=paste0(prefix, "_auto.png"), height=4, width=8, units="in", dpi=300)

  p1 <- make_sim_plot(tot, names=c("fretHuberauto300", "fretPoisauto300", "fretTauto300"),
                      lty=c(1, 2, 3), col=rep("black", 3), shapes=c(15, 16, 17))
  h <- grid.arrange(p1$fdpplot, p1$tprplot, ncol=2)
  ggsave(h, file=paste0(prefix, "_fret.png"), height=4, width=8, units="in", dpi=300)

  p1 <- make_sim_plot(tot, names=c("fretHuberauto150", paste0(c("Pois", "Huber", "T"), "naive64")),
                      col=c("black", "violetRed", "chartreuse3", "blue"),
                      lty=c(1, 2, 1,3), shapes=15:18)
  h <- grid.arrange(p1$fdpplot, p1$tprplot, ncol=2)
  ggsave(h, file=paste0(prefix, "_fretauto150_naive64.png"), height=4, width=8, units="in", dpi=300)

  p1 <- make_sim_plot(tot, names=c("fretHuberauto50", paste0(c("Pois", "Huber", "T"), "informed32")),
                      col=c("black", "violetRed", "chartreuse3", "blue"),
                      lty=c(1, 2, 1,3), shapes=15:18)
  h <- grid.arrange(p1$fdpplot, p1$tprplot, ncol=2)
  ggsave(h, file=paste0(prefix, "_fretauto50_informed32.png"), height=4, width=8, units="in", dpi=300)


  #l <- make_sim_legend_aoas()
  #ggsave(l$plot, file="~/Dropbox/Thesis/final_talk/img/legend.png",
  #       height=4, width=3, units="in", dpi=300)



}


make_sim_legend_aoas <- function(){
  points <- data.frame(x=rep(1, 5), y = rev(seq(1, 3, length.out=5)),
                       x2 = 1:5, y2=rep(1, 5))

  points$left = points$x2-0.45
  points$right = points$x2 + 0.45

  points$labs <- c("FRET", "DESeq2", "WaveQTL", "Huber", "QP")
  points$lty <- c(1, 2, 4, 2, 3)
  points$shape <- 15:19
  points$informed <- c("fretHuber60", paste0("w64b", c("DESeq2", "Wave", "Huber", "Pois")))
  points$naive <- c("fretHuber60", paste0("w64e", c("DESeq2", "Wave", "Huber", "Pois")))

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

