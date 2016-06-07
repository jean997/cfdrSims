
merge_fig <- function(){
  D <- getobj("~/Desktop/Cluster_FDR/region_51/huber_stats_r51.RData")
x <- D[,2]
pos <- D[,1]
rm(D)
xs <- ksmooth(x=pos, y=x, bandwidth=20, x.points = pos)$y
z0 = quantile(abs(xs), probs=0.9)*0.3
xxs <- xs[50400:52000]
p <- pos[50400:52000]
dat <- data.frame(cbind(p, xxs))
n <- dim(dat)[1]
dat$col <- rep("black", n)
dat$col[abs(xxs) > 1.1] <- "blue"

png("~/Dropbox/CFDR/img/merge.png", height=4, width=5, units = "in", res=300)
plot(dat$p, dat$xxs, type="n", xlab="Position", ylab="T(j)")
segments(x0=dat$p[-n], y0=dat$xxs[-n], x1=dat$p[-1], y1=dat$xxs[-1], col=dat$col, lty=1)
abline(h=-0.4, lty=2, lwd=1.5)
segments(x0=c(dat$p[1], dat$p[550], dat$p[802]), y0=rep(-1.1, 3),
         x1=c(dat$p[549], dat$p[801], dat$p[n]), y1=rep(-1.1, 3),
         col=c("black", "blue", "black"), lwd=1.5)

dev.off()
}
