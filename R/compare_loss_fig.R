

compare_loss_fig <- function(c=1, lim=c(-3, 3), npoints=1000){
  x <- seq(lim[1], lim[2], length.out=npoints)
  y1 <- 0.5*x^2
  y2 <- abs(x)
  y3 <- y1
  y3[ abs(x) > c] <- c*(abs(x[abs(x) > c]) - 0.5*c)
  d <- data.frame(cbind(x, y1, y2, y3))
  names(d) <- c("x", "Squared Error", "Absolute Value", paste0("Huber c=", c))
  dlong <- gather(d, Loss, y, -x)
  h <- ggplot(dlong) + geom_line(aes(x=x, y=y, group=Loss, lty=Loss, col=Loss), lwd=1.2) +
    labs(x="u", y=bquote(rho~"(u)")) + ylim(c(0, lim[2])) +
    theme_bw() + theme(panel.grid=element_blank())
  return(h)
}
