

make_plot <- function(summ.files, out.name, lty.legend, ltys=c(2, 1), ymax=NULL){

  p = length(summ.files)

  stopifnot(length(ltys) == p)

  cols=c("seagreen", "violetRed", "black")

  png(out.name, width=1000, height=500)
  par(mfrow=c(1, 2))
  if(is.null(ymax)){
    ymax <- 0.2
    for(file in summ.files){
      Z <- getobj(file)
      fdp <- apply(Z$fdp, MARGIN=c(1, 2), FUN=mean)
      ymax <- max(ymax, max(fdp))
    }
  }
  plot(0, 0, yaxt="n", xaxt="n", type="n", xlab="Target FDR",
           ylab="Average False Discovery Proportion", xlim=c(0, 0.2),
           ylim=c(0, ymax))
  axis(side = 1, at=c(0.02, 0.05, 0.1, 0.2))
  axis(side = 2, at=c(0.02, 0.05, 0.1, 0.2))
  abline(0, 1, lty=2)

  for(j in 1:p){
    Z <- getobj(summ.files[j])
    fdp <- apply(Z$fdp, MARGIN=c(1, 2), FUN=mean)
    for(i in 1:3){
      lines(c(0.02, 0.05, 0.1, 0.2), fdp[i,], lty=ltys[j],
            col=cols[i], pch=i, type="b", lwd=3)
    }
  }

  plot(0, 0, yaxt="n", xaxt="n", type="n", xlab="Target FDR",
       ylab="Average True Positive Rate", xlim=c(0, 0.2),
       ylim=c(0,1))
  axis(side = 1, at=c(0.02, 0.05, 0.1, 0.2))
  axis(side = 2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1))

  for(j in 1:p){
    Z <- getobj(summ.files[j])
    tpr <- apply(Z$tpr, MARGIN=c(1, 2), FUN=mean)
    for(i in 1:3){
      lines(c(0.02, 0.05, 0.1, 0.2), tpr[i,], lty=ltys[j],
            col=cols[i], pch=i, type="b", lwd=3)
    }

  }
  legend("topleft", legend=c("Poisson", "Huber", "t-test"),
         col=c("black", "violetRed", "seagreen"), lty=1)
  #if(length(s0s) ==2){
  #  l <- c()
  #  for(j in 1:p){
  #    l <- c(l, as.expression(bquote(delta == .(s0s[j]))))
  #  }
  #  legend("topright", legend=l, lty=c(2, 1))
  #}else{
  #  l <- c()
  #  for(j in 1:p){
  #    l <- c(l, as.expression(bquote("nsegs" == .(nsegs[j]))))
  #  }
  legend("topright", legend=lty.legend, lty=ltys)

  dev.off()

}
