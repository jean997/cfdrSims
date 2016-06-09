

plot_performance <- function(){
  R <- getobj("sim_results/ss916_s0_0.05_window_rates.RData")
  cfdr05 <- avg_by_ct(R$tpp[[1]], R$fp[[1]])
  cfdr05_fdp <- avg_by_fdp(R$tpp[[1]], R$fdp[[1]])

  R <- getobj("sim_results/ss916_s0_0_window_rates.RData")
  fdp <- apply(R$fdp_at, MARGIN=c(2, 3), FUN=mean)
  tpr <- apply(R$tpr_at, MARGIN=c(2, 3), FUN=mean)

  cfdr0 <- avg_by_ct(R$tpp[[1]], R$fp[[1]])
  opt20 <- avg_by_ct(R$tpp[[2]], R$fp[[2]])
  bin50 <- avg_by_ct(R$tpp[[4]], R$fp[[4]])
  bin20 <- avg_by_ct(R$tpp[[6]], R$fp[[6]])

  cfdr0_fdp <- avg_by_fdp(R$tpp[[1]], R$fdp[[1]])
  opt20_fdp <- avg_by_fdp(R$tpp[[2]], R$fdp[[2]])
  bin50_fdp <- avg_by_fdp(R$tpp[[4]], R$fdp[[4]])
  bin20_fdp <- avg_by_fdp(R$tpp[[6]], R$fdp[[6]])

  R <- getobj("sim_results/ss916_s0_0_window_rates_2.RData")
  fdp_2 <- apply(R$fdp_at, MARGIN=c(2, 3), FUN=mean)
  tpr_2 <- apply(R$tpr_at, MARGIN=c(2, 3), FUN=mean)

  opt30 <- avg_by_ct(R$tpp[[1]], R$fp[[1]])
  opt50 <- avg_by_ct(R$tpp[[2]], R$fp[[2]])
  opt100 <- avg_by_ct(R$tpp[[3]], R$fp[[3]])
  regions <- avg_by_ct(R$tpp[[4]], R$fp[[4]])
  opt30_fdp <- avg_by_fdp(R$tpp[[1]], R$fdp[[1]])
  opt50_fdp <- avg_by_fdp(R$tpp[[2]], R$fdp[[2]])
  opt100_fdp <- avg_by_fdp(R$tpp[[3]], R$fdp[[3]])
  regions_fdp <- avg_by_fdp(R$tpp[[4]], R$fdp[[4]])

  R <- getobj("sim_results/ss916_s0_0_waveQTL.RData")
  fdp_w <- apply(R$fdp_at, MARGIN=c(2, 3), FUN=mean)
  tpr_w <- apply(R$tpr_at, MARGIN=c(2, 3), FUN=mean)
  wave64 <- avg_by_ct(R$tpp[[1]], R$fp[[1]])
  wave128 <- avg_by_ct(R$tpp[[2]], R$fp[[2]])
  wave64_fdp <- avg_by_fdp(R$tpp[[1]], R$fdp[[1]])
  wave128_fdp <- avg_by_fdp(R$tpp[[2]], R$fdp[[2]])

  #Plot tpr vs. number of false postives
  plot(cfdr05$fct, cfdr05$tpr, ylim=c(0, 1), xlim=c(0, 40), type="l", xlab="False Positives", ylab="True Positive Rate")
  lines(cfdr0$fct, cfdr0$tpr, lty=2)
  lines(opt20$fct, opt20$tpr, col="blue")
  #lines(opt30$fct, opt30$tpr, col="blue", lty=2)
  lines(opt50$fct, opt50$tpr, col="orange")
  #lines(opt100$fct, opt100$tpr, col="orange", lty=2)
  lines(regions$fct, regions$tpr, col="red")
  lines(bin20$fct, bin20$tpr, col="purple")
  lines(bin50$fct, bin50$tpr, col="purple", lty=2)
  lines(wave64$fct, wave64$tpr, col="seagreen")
  lines(wave128$fct, wave128$tpr, col="seagreen", lty=2)

  legend("bottomright", legend=c("discopony", "opt", "opt50", "opt200", "bin20", "bin50", "wave64", "wave128"),
         col=c("black", "blue", "orange", "red", "purple", "purple", "seagreen", "seagreen"),
         lty=c(1, 1, 1, 1, 2, 1, 2))

  #Plot max tpr vs. false discovery rate
  plot(cfdr05_fdp$fdp, cfdr05_fdp$tpr, ylim=c(0, 1), xlim=c(0, 1), type="l",
       xlab="False Discovery Rate", ylab="Max True Positive Rate")
  lines(cfdr0_fdp$fdp, cfdr0_fdp$tpr, lty=2)
  lines(opt20_fdp$fdp, opt20_fdp$tpr, col="blue")
  lines(opt50_fdp$fdp, opt50_fdp$tpr, col="orange")
  lines(regions_fdp$fdp, regions_fdp$tpr, col="red")
  lines(bin20_fdp$fdp, bin20_fdp$tpr, col="purple")
  lines(bin50_fdp$fdp, bin50_fdp$tpr, col="purple", lty=2)
  lines(wave64_fdp$fdp, wave64_fdp$tpr, col="seagreen")
  lines(wave128_fdp$fdp, wave128_fdp$tpr, col="seagreen", lty=2)
  abline(v=0.1)
  legend("bottomright", legend=c("discopony", "opt", "opt50", "opt200", "bin20", "bin50", "wave64", "wave128"),
         col=c("black", "blue", "orange", "red", "purple", "purple", "seagreen", "seagreen"),
         lty=c(1, 1, 1, 1, 2, 1, 2))


  #Plot fdr at nominal levels
  levels=c(0.02, 0.05, 0.1, 0.2)
  plot(0, 0, yaxt="n", xaxt="n", type="n", xlab="Target FDR",
       ylab="Average False Discovery Proportion", xlim=c(0, max(levels)),
       ylim=c(0, 0.45))
  axis(side = 1, at=levels)
  axis(side = 2, at=levels)
  abline(0, 1, lty=2)
  lines(levels, fdp[1,],  pch=1, type="b", lwd=1.5, col="black") #discopony
  lines(levels, fdp[2,],  pch=2, type="b", lwd=1.5, col="blue") #opt20
  lines(levels, fdp[6,],  pch=3, type="b", lwd=1.5, col="purple") #bin20
  lines(levels, fdp[4,],  pch=3, type="b", lwd=1.5, col="purple", lty=2) #bin50
  lines(levels, fdp_2[2,],  pch=2, type="b", lwd=1.5, col="orange") #Opt50
  lines(levels, fdp_2[4,],  pch=2, type="b", lwd=1.5, col="red") #regions
  lines(levels, fdp_w[1,],  pch=4, type="b", lwd=1.5, col="seagreen") #wave64
  lines(levels, fdp_w[2,],  pch=4, type="b", lwd=1.5, col="seagreen", lty=2) #wave64

  #Plot tpr at nominal levels
  plot(0, 0, yaxt="n", xaxt="n", type="n", xlab="Target FDR",
       ylab="Average True Positive Rate", xlim=c(0, max(levels)),
       ylim=c(0,1))
  axis(side = 1, at=levels)
  axis(side = 2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1))
  lines(levels, tpr[1,],  pch=1, type="b", lwd=1.5, col="black") #discopony
  lines(levels, tpr[2,],  pch=2, type="b", lwd=1.5, col="blue") #opt20
  lines(levels, tpr[6,],  pch=3, type="b", lwd=1.5, col="purple") #bin20
  lines(levels, tpr[4,],  pch=3, type="b", lwd=1.5, col="purple", lty=2) #bin50
  lines(levels, tpr_2[2,],  pch=2, type="b", lwd=1.5, col="orange") #Opt50
  lines(levels, tpr_2[4,],  pch=2, type="b", lwd=1.5, col="red") #regions
  lines(levels, tpr_w[1,],  pch=4, type="b", lwd=1.5, col="seagreen") #wave64
  lines(levels, tpr_w[2,],  pch=4, type="b", lwd=1.5, col="seagreen", lty=2) #wave64



}
