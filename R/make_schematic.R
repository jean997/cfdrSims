make_schematic <- function(){

  w=30
  ss <- seq((w/20), w, length.out=10)
  iv <- w/20
  cex=1
plot(0, 0, type="n", ylab="", xlab="", xaxt="n",
     yaxt="n", bty="n", ylim=c(0, 0.5), xlim=c(0, 6*(w+iv)+1))
rect(0, 0.2, w + iv, 0.3, lwd=2)
text(x=ss, y=rep(0.25, 10), labels = "1", cex=cex)
rect(w + iv, 0.2, 2*(w+iv), 0.3, lwd=2)
text(x=ss+w+iv, y=rep(0.25, 10), labels = "*", cex=cex)
rect(2*(w+iv), 0.2, 3*(w+iv), 0.3, lwd=2)
text(x=ss +(2*(w+iv)), y=rep(0.25, 10), labels = "1", cex=cex)

rect(3*(w+iv), 0.2, 4*(w+iv), 0.3, lwd=2)
text(x=ss+(3*(w+iv)), y=rep(0.25, 10), labels = "2", cex=cex)
rect(4*(w+iv), 0.2, 5*(w+iv), 0.3, lwd=2)
text(x=ss+(4*(w+iv)), y=rep(0.25, 10), labels = "*", cex=cex)
rect(5*(w+iv), 0.2, 6*(w+iv), 0.3, lwd=2)
text(x=ss+(5*(w+iv)), y=rep(0.25, 10), labels = "2", cex=cex)

}


make_schematic_talk <- function(){

  w=30
  ss <- seq((w/10), w, length.out=10)

  rect.data <- data.frame(xmin = (0:5)*w, xmax=(1:6)*w,
                          ymin=rep(0, 6), ymax=rep(0.3, 6),
                          type=factor(c(1, 3, 1, 2, 3, 2)))
  seg.data <- data.frame(x=seq(w/10, 6*w, length.out=60), xend=seq(w/10, 6*w, length.out=60),
                        y = rep(0, 60), yend=rep(0.3, 60))
  h <- ggplot(rect.data) + geom_rect(aes(xmin=xmin, xmax=xmax,
                                         ymin=ymin, ymax=ymax, fill=type), lwd=1, col="black") +
    geom_segment(data=seg.data, aes(x=x, xend=xend, y=y, yend=yend)) +
    theme_bw(18) + theme( panel.grid=element_blank(), panel.border = element_blank(),
                           axis.text=element_blank(), axis.ticks=element_blank(),
                             legend.position="none", axis.title=element_blank())

  ggsave("~/Dropbox/Thesis/final_talk/img/schematic.png", height=1, width=7, units="in", dpi=300)

  plot(0, 0, type="n", ylab="", xlab="", xaxt="n",
       yaxt="n", bty="n", ylim=c(0, 0.5), xlim=c(0, 6*(w+iv)+1))
  rect(0, 0.2, w + iv, 0.3, lwd=2, col="blue")
  abline(v=ss)
  text(x=ss, y=rep(0.25, 10), labels = "1", cex=cex)
  rect(w + iv, 0.2, 2*(w+iv), 0.3, lwd=2)
  text(x=ss+w+iv, y=rep(0.25, 10), labels = "*", cex=cex)
  rect(2*(w+iv), 0.2, 3*(w+iv), 0.3, lwd=2)
  text(x=ss +(2*(w+iv)), y=rep(0.25, 10), labels = "1", cex=cex)

  rect(3*(w+iv), 0.2, 4*(w+iv), 0.3, lwd=2)
  text(x=ss+(3*(w+iv)), y=rep(0.25, 10), labels = "2", cex=cex)
  rect(4*(w+iv), 0.2, 5*(w+iv), 0.3, lwd=2)
  text(x=ss+(4*(w+iv)), y=rep(0.25, 10), labels = "*", cex=cex)
  rect(5*(w+iv), 0.2, 6*(w+iv), 0.3, lwd=2)
  text(x=ss+(5*(w+iv)), y=rep(0.25, 10), labels = "2", cex=cex)

}
