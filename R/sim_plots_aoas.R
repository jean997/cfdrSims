
sim_plots_aoas <- function(){
  #tot.rates <- getobj("sim_results/bin1_total3.RData")
  tot.rates <- getobj("sim_results/bin3_64.RData")
  tot.rates <- getobj("sim_results/bin1_64.RData")
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

  p1 <- make_sim_plot(tot.rates, names=c("fretHuber2", "fretHuber1", "fretHuber20", "fretHuber60"),
                      col=c("black", "violetRed", "chartreuse3", "blue"),
                      lty=rep(1, 4), shapes=rep(15, 4))
  #ggsave(p1$tprplot, file="~/Dropbox/CFDR-Jean/for_AOAS/img/fret_compare2_tpp.png",height=4, width=4, units="in", dpi=300)
  #ggsave(p1$fdpplot, file="~/Dropbox/Thesis/final_talk/img/fret_compare2_fdp.png",height=4, width=4, units="in", dpi=300)

  p2 <- make_sim_plot(tot.rates, names=c("fretHuber60", "fretPois60", "fretT60"),
                      lty=c(1, 2, 3), col=rep("black", 3), shapes=c(15, 16, 17))
  #ggsave(p2$tprplot, file="~/Dropbox/Thesis/final_talk/img/fret_compare_tpp.png",height=4, width=4, units="in", dpi=300)
  #ggsave(p2$fdpplot, file="~/Dropbox/Thesis/final_talk/img/fret_compare_fdp.png",height=4, width=4, units="in", dpi=300)

  l <- make_sim_legend_aoas()
  #ggsave(l$plot, file="~/Dropbox/Thesis/final_talk/img/legend.png",
  #       height=4, width=3, units="in", dpi=300)

  p64b <- make_sim_plot(tot.rates, names=l$info$informed,
                        cols=l$info$cols,
                        ltys=l$info$lty, shapes=l$info$shape)
  #ggsave(p64b$fdpplot, file="~/Dropbox/Thesis/final_talk/img/w64b_fdp.png", height=4, width=4, units="in", dpi=300)
  #ggsave(p64b$tprplot, file="~/Dropbox/Thesis/final_talk/img/w64b_tpp.png", height=4, width=4, units="in", dpi=300)


  p64e <- make_sim_plot(tot.rates, names=l$info$naive,
                        cols=l$info$cols,
                        ltys=l$info$lty, shapes=l$info$shape)
  #ggsave(p64e$fdpplot, file="~/Dropbox/Thesis/final_talk/img/w64e_fdp.png", height=4, width=4, units="in", dpi=300)
  #ggsave(p64e$tprplot, file="~/Dropbox/Thesis/final_talk/img/w64e_tpp.png", height=4, width=4, units="in", dpi=300)


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

