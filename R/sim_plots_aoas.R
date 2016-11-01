
sim_plots_aoas <- function(){
  tot <- getobj("sim_results/bin5_pb64_total.RData")

  p1 <- make_sim_plot(tot, names=c("fretHuberauto3200",
                                   "fretPoisauto3200",
                                   "DESeq2informed64",
                                   "Waveinformed64",
                                   "Huberinformed64"),
                      col=c("black", "black",
                            "violetRed", "chartreuse3", "blue"),
                      lty=c(1, 2, 3, 4, 2), shapes=c(16, 15, 17, 18, 3))
  h <- grid.arrange(p1$fdpplot, p1$tprplot, ncol=2)
  ggsave(h, file="~/Dropbox/cfdr-jean/for_AOAS/img/bin5_informed.png", height=4, width=8, units="in", dpi=300)

  p2 <- make_sim_plot(tot, names=c("fretHuberauto3200",
                                   "fretPoisauto3200",
                                   "DESeq2naive64",
                                   "Wavenaive64",
                                   "Hubernaive64"),
                      col=c("black", "black","violetRed", "chartreuse3", "blue"),
                      lty=c(1, 2, 3, 4, 2), shapes=c(16, 15, 17, 18, 3))
  h <- grid.arrange(p2$fdpplot, p2$tprplot, ncol=2)
  ggsave(h, file="~/Dropbox/cfdr-jean/for_AOAS/img/bin5_naive.png", height=4, width=8, units="in", dpi=300)

  h <- grid.arrange(p1$fdpplot + ggtitle("FDR; Informed Bins"), p1$tprplot+ggtitle("Power; Informed Bins"),
                    p2$fdpplot + ggtitle("FDR; Naive Bins"), p2$tprplot+ggtitle("Power; Naive Bins"))
  ggsave(h, file="~/Dropbox/cfdr-jean/for_AOAS/img/bin5.png", height=7, width=8, units="in", dpi=300)


  nn <- c("DESeq2informed64", "DESeq2naive64", "Huberinformed64", "Hubernaive64",
          "Waveinformed64", "Wavenaive64")
  ttls <- c("DESeq2 Informed Bins", "DESeq2 Naive Bins",
            "Huber Informed Bins", "Huber Naive Bins",
            "WaveQTL Informed Bins", "WaveQTL Naive Bins")
  h <- list(); ct <- 1
  for(i in 1:length(nn)){
    ix <- which(tot$names==nn[i])
    null.p <- tot$stats[[ix]]$pval[tot$stats[[ix]]$signal==0]
    h[[ct]] <- thin_qqplot(null.p)  + ggtitle(ttls[i])
    ct <- ct + 1
  }
  hh <- grid.arrange(h[[1]], h[[2]], h[[3]], h[[4]], h[[5]], h[[6]], ncol=2)
  ggsave(hh, file="~/Dropbox/cfdr-jean/for_AOAS/img/sim_qqplots.png", height=12, width=8, units="in", dpi=300)

  #width 32 bins
  p1 <- make_sim_plot(tot, names=c("fretHuberauto3200",
                                   "fretPoisauto3200",
                                   "DESeq2informed32",
                                   "Waveinformed32",
                                   "Huberinformed32"),
                      col=c("black", "black",
                            "violetRed", "chartreuse3", "blue"),
                      lty=c(1, 2, 3, 4, 2), shapes=c(16, 15, 17, 18, 3))

  p2 <- make_sim_plot(tot, names=c("fretHuberauto3200",
                                   "fretPoisauto3200",
                                   "DESeq2naive32",
                                   "Wavenaive32",
                                   "Hubernaive32"),
                      col=c("black", "black","violetRed", "chartreuse3", "blue"),
                      lty=c(1, 2, 3, 4, 2), shapes=c(16, 15, 17, 18, 3))

  h <- grid.arrange(p1$fdpplot + ggtitle("FDR; Informed Bins"), p1$tprplot+ggtitle("Power; Informed Bins"),
                    p2$fdpplot + ggtitle("FDR; Naive Bins"), p2$tprplot+ggtitle("Power; Naive Bins"))
  ggsave(h, file="~/Dropbox/cfdr-jean/for_AOAS/img/bin5_32.png", height=7, width=8, units="in", dpi=300)

  h <- list(); ct <- 1
  nn <- c("DESeq2informed32", "DESeq2naive32", "Huberinformed32", "Hubernaive32",
          "Waveinformed32", "Wavenaive32")
  for(i in 1:length(nn)){
    ix <- which(tot$names==nn[i])
    null.p <- tot$stats[[ix]]$pval[tot$stats[[ix]]$signal==0]
    h[[ct]] <- thin_qqplot(null.p)  + ggtitle(ttls[i])
    ct <- ct + 1
  }
  hh <- grid.arrange(h[[1]], h[[2]], h[[3]], h[[4]], h[[5]], h[[6]], ncol=2)
  ggsave(hh, file="~/Dropbox/cfdr-jean/for_AOAS/img/sim_qqplots_32.png", height=12, width=8, units="in", dpi=300)



  #Width 128 bins
  p1 <- make_sim_plot(tot, names=c("fretHuberauto3200",
                                   "fretPoisauto3200",
                                   "DESeq2informed128",
                                   "Waveinformed128",
                                   "Huberinformed128"),
                      col=c("black", "black",
                            "violetRed", "chartreuse3", "blue"),
                      lty=c(1, 2, 3, 4, 2), shapes=c(16, 15, 17, 18, 3))

  p2 <- make_sim_plot(tot, names=c("fretHuberauto3200",
                                   "fretPoisauto3200",
                                   "DESeq2naive128",
                                   "Wavenaive128",
                                   "Hubernaive128"),
                      col=c("black", "black","violetRed", "chartreuse3", "blue"),
                      lty=c(1, 2, 3, 4, 2), shapes=c(16, 15, 17, 18, 3))

  h <- grid.arrange(p1$fdpplot + ggtitle("FDR; Informed Bins"), p1$tprplot+ggtitle("Power; Informed Bins"),
                    p2$fdpplot + ggtitle("FDR; Naive Bins"), p2$tprplot+ggtitle("Power; Naive Bins"))
  ggsave(h, file="~/Dropbox/cfdr-jean/for_AOAS/img/bin5_128.png", height=7, width=8, units="in", dpi=300)

  h <- list(); ct <- 1
  nn <- c("DESeq2informed128", "DESeq2naive128", "Huberinformed128", "Hubernaive128",
          "Waveinformed128", "Wavenaive128")
  for(i in 1:length(nn)){
    ix <- which(tot$names==nn[i])
    null.p <- tot$stats[[ix]]$pval[tot$stats[[ix]]$signal==0]
    h[[ct]] <- thin_qqplot(null.p)  + ggtitle(ttls[i])
    ct <- ct + 1
  }
  hh <- grid.arrange(h[[1]], h[[2]], h[[3]], h[[4]], h[[5]], h[[6]], ncol=2)
  ggsave(hh, file="~/Dropbox/cfdr-jean/for_AOAS/img/sim_qqplots_128.png", height=12, width=8, units="in", dpi=300)



}


make_sim_legend_aoas <- function(){
  #points <- data.frame(x=rep(1, 5), y = rev(seq(1, 3, length.out=5)))
  points <- data.frame(y=rep(1, 5), x = seq(1, 2.25, length.out=5))

  points$left = points$x-0.12
  points$right = points$x + 0.12

  points$labs <- c("FRET-Huber", "FRET-QP", "DESeq2", "WaveQTL", "Huber")
  points$lty <-c(1, 2, 3, 4, 2)
  points$shape <- c(16, 15, 17, 18, 3)

  points$cols <- c("black", "black","violetRed", "chartreuse3", "blue")

  p <-ggplot(points) +
    geom_segment(aes(x=left, xend=right, y=y, yend=y), lty=points$lty,
                 colour=points$cols, lwd=1.3)+
    geom_point(aes(x=x, y=y), col="white", shape=20, size=13)+
    geom_point(aes(x=x, y=y), col=points$cols, shape=points$shape, size=2, stroke=1.3)+
    annotate(geom="text", y=rep(1.1, 5), x=points$x,
             label=points$labs, size=4) +
    xlim(0.8, 2.5) + ylim(0.97, 1.2) +
    theme_bw() + theme(panel.grid=element_blank(), axis.title=element_blank(),
                       axis.text=element_blank(),
                       panel.border=element_blank(), axis.ticks=element_blank())

  ggsave(p, file="~/Dropbox/cfdr-jean/for_AOAS/img/sim_legend.png", height=1, width=8, units="in", dpi=300)
  return(list("plot"=p, "info"=points))
}

example_profile_aoas <- function(){
  set.seed(1e8)
  peak.starts <- sort(sample(1:800, size=4))

  peak.hts <- rexp(n=4, 1/5)
  yy <- cfdrSims:::gen_profile2(peak.starts, peak.hts, total.length = 800,
                                peak.base=64, mesa.width=64/3,
                                bw=64/4, bg.ht=1.5)

  df <- data.frame(pos=1:800, y=yy)

  p <- 800
  k <- 4

  informed <- data.frame(cbind(peak.starts, peak.starts + 63))
  names(informed) = c("x", "xend")
  informed$type <- "informed"
  informed$y <- 1

  naive <- data.frame(cbind(seq(1, p - (p%%64), by=64), seq(64, p, by=64)))
  names(naive)=c("x", "xend")
  naive$type="naive"
  naive$y <-rep(c(0.5, 0.6), 6)
  allbins <- rbind(informed, naive)

  h <- ggplot(df) + geom_line(aes(x=pos, y=yy)) +
    geom_segment(data=allbins, aes(x=x, xend=xend, y=y, yend=y, col=type), lwd=1.5) +
    #scale_color_manual(labels=c("Informed Bins", "Naive Bins"), values=c("blue", "darkorange2")) +
    #scale_y_continuous(breaks=c(1.5, 2:7)) + ylab(expression("Profile ("~gamma~")")) +
    xlab("Position") + ylab(expression("Profile ("~gamma~")")) +
    theme_bw(12) + theme(#axis.title.x=element_blank(),
      panel.grid=element_blank(),
      legend.position="none")
  ggsave(h, file="~/Dropbox/cfdr-jean/for_AOAS/img/example_profile.png", height=2.5, width=8, units="in", dpi=300)

}
