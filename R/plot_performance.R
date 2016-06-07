

plot_performance <- function(){
  load("sim_results/ss916_s0_0.05_window_rates.RData")
  cfdr05 <- avg_by_ct(tpr.list=tpp[[1]], fct.list = fp[[1]])
  cfdr05_fdp <- avg_by_interp(tpp[[1]], fdp[[1]])

  load("sim_results/ss916_s0_0_window_rates.RData")
  cfdr0 <- avg_by_ct(tpr.list=tpp[[1]], fct.list = fp[[1]])
  opt20 <- avg_by_ct(tpp[[2]], fp[[2]])
  bin50 <- avg_by_ct(tpp[[4]], fp[[4]])
  bin20 <- avg_by_ct(tpp[[6]], fp[[6]])

  cfdr0_fdp <- avg_by_interp(tpp[[1]], fdp[[1]])
  opt20_fdp <- avg_by_interp(tpp[[2]], fdp[[2]])
  bin50_fdp <- avg_by_interp(tpp[[4]], fdp[[4]])
  bin20_fdp <- avg_by_interp(tpp[[6]], fdp[[6]])

  load("sim_results/ss916_s0_0_window_rates_2.RData")
  opt30 <- avg_by_ct(tpp[[1]], fp[[1]])
  opt50 <- avg_by_ct(tpp[[2]], fp[[2]])
  opt100 <- avg_by_ct(tpp[[3]], fp[[3]])
  regions <- avg_by_ct(tpp[[4]], fp[[4]])
  opt30_fdp <- avg_by_interp(tpp[[1]], fdp[[1]])
  opt50_fdp <- avg_by_interp(tpp[[2]], fdp[[2]])
  opt100_fdp <- avg_by_interp(tpp[[3]], fdp[[3]])
  regions_fdp <- avg_by_interp(tpp[[4]], fdp[[4]])


  plot(cfdr05$fct, cfdr05$tpr, ylim=c(0, 1), xlim=c(0, 40), type="l", xlab="False Positives", ylab="True Positive Rate")
  lines(cfdr0$fct, cfdr0$tpr, lty=2)
  lines(opt20$fct, opt20$tpr, col="blue")
  #lines(opt30$fct, opt30$tpr, col="blue", lty=2)
  lines(opt50$fct, opt50$tpr, col="orange")
  #lines(opt100$fct, opt100$tpr, col="orange", lty=2)
  lines(regions$fct, regions$tpr, col="red")
  lines(bin20$fct, bin20$tpr, col="purple")
  lines(bin50$fct, bin50$tpr, col="purple", lty=2)

  legend("bottomright", legend=c("discopony", "opt", "opt50", "opt200", "bin20", "bin50"),
         col=c("black", "blue", "orange", "red", "purple", "purple"),
         lty=c(1, 1, 1, 1, 2))

  plot(cfdr05_fdp$fpr, cfdr05_fdp$tpr, ylim=c(0, 1), xlim=c(0, 1), type="l",
       xlab="False Discovery Proportion", ylab="True Positive Rate")
  lines(cfdr0_fdp$fpr, cfdr0_fdp$tpr, lty=2)
  lines(opt20_fdp$fpr, opt20_fdp$tpr, col="blue")
  lines(opt50_fdp$fpr, opt50_fdp$tpr, col="orange")
  lines(regions$fpr, regions$tpr, col="red")
  lines(bin20_fdp$fpr, bin20_fdp$tpr, col="purple")
  lines(bin50_fdp$fpr, bin50_fdp$tpr, col="purple", lty=2)
  abline(v=0.1)
  legend("bottomright", legend=c("discopony", "opt", "opt50", "opt200", "bin20", "bin50"),
         col=c("black", "blue", "orange", "red", "purple", "purple"),
         lty=c(1, 1, 1, 1, 2))

}
