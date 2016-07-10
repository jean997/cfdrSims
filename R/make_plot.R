
#' Plot Simulations
#'@description Plot simulations
#'@param summ.files List of summary files (see collect.R)
#'@param lty.legend Legend labels for line type
#'@param ltys Line type for each summary file
#'@param ymax Optional max height of FDR y axis
#'@param levels Levels in summ.files
#' @return nothing
#'@export
make_sim_plot <- function(tot.rates, names, ltys, cols, shapes,
                      levels=c(0.02, 0.05, 0.1, 0.2)){


  p = length(names)
  stopifnot(length(ltys) == p)
  if(is.null(cols)) cols=1:p

  stopifnot(all(names %in% tot.rates$names))
  ix <- match(names, tot.rates$names)

  tpr <- tot.rates$avg.tpr[ix,]
  tpr <- data.frame(tpr)
  fdp <- tot.rates$avg.fdp[ix,]
  fdp <- data.frame(fdp)

  tpr$type <- fdp$type <- names
  names(tpr) <- names(fdp) <- c(levels, "type")
  tprlong <- gather(tpr, "level", "tpr", -type, convert=TRUE)
  fdplong <- gather(fdp, "level", "fdp", -type, convert=TRUE)
  tprlong$type = factor(tprlong$type, levels = names)
  fdplong$type = factor(fdplong$type, levels = names)
  #tprlong$type <- relevel(tprlong$type, names)
  #fdplong$type <- relevel(fdplong$type, names)

  tprplot <- ggplot(tprlong) +
    geom_line(aes(x=level, y=tpr, group=type, col=type, lty=type), lwd=1.3) +
    geom_point(aes(x=level, y=tpr), size=5, colour="white", shape=20) +
    geom_point(aes(x=level, y=tpr, col=type, shape=type), size=2, stroke=1.3) +
    theme_bw() + labs(x="Target FDR", y="Average True Positive Rate") +
    scale_color_manual(values=cols) +
    scale_linetype_manual(values=ltys) +
    scale_shape_manual(values=shapes) + ggtitle("True Positive Rate") +
    scale_x_continuous(breaks=levels, limits=c(0, max(levels))) +
    scale_y_continuous(limits = c(0, max(tot.rates$avg.tpr[ix,]))) +
    theme(panel.grid=element_blank(), legend.position="none")

  if(max(tot.rates$avg.fdp[ix,]) < max(levels)) ymx <- max(levels)
    else ymx <- max(tot.rates$avg.fdp[ix,])
  fdpplot <- ggplot(fdplong) +
    geom_line(aes(x=level, y=fdp, group=type, col=type, lty=type), lwd=1.3) +
    geom_point(aes(x=level, y=fdp), size=5, colour="white", shape=20) +
    geom_point(aes(x=level, y=fdp, col=type, shape=type), size=2, stroke=1.3) +
    geom_abline(slope=1, intercept=0) +
    theme_bw() + labs(x="Target FDR", y="Average False Discovery Proportion") +
    scale_color_manual(values=cols) +
    scale_linetype_manual(values=ltys) +
    scale_shape_manual(values=shapes) + ggtitle("False Discovery Rate") +
    scale_x_continuous(breaks=levels, limits=c(0, max(levels))) +
    scale_y_continuous(breaks=levels, limits = c(0, ymx)) +
    theme(panel.grid=element_blank(), legend.position="none")
  return(list("tprplot"=tprplot, "fdpplot"=fdpplot))
}

make_sim_legend <- function(){
  points <- cbind(rep(c(1, 2, 3), each=3), rep(c(1, 1.5, 2.3), 3))
  points <- data.frame(points)
  names(points)=c("x", "y")
  points$left = points$x-0.45
  points$right = points$x + 0.45

  points$lty = 1
  points$lty[points$y==1] <- 2
  points$lty[points$x==2 & points$y==2.3] <- 2
  points$lty[points$x==3 & points$y==2.3] <- 3

  points$cols <- "black"
  points$cols[points$x==1 & points$y < 2] <- "blue"
  points$cols[points$x==2 & points$y < 2] <- "violetRed"
  points$cols[points$x==3 & points$y < 2] <- "chartreuse3"

  points$shape = c(1, 1, 15, 2, 2, 16, 3, 3, 17)

  points$name <- NA
  points$name[points$col=="black"] <- paste0("fretHuber", c(1, 2, 6))
  points$name[points$cols=="blue"] <- paste0(c("w64e", "w64b"), "DESeq2")
  points$name[points$cols=="violetRed"] <- paste0(c("w64e", "w64b"), "Wave")
  points$name[points$cols=="chartreuse3"] <- paste0(c("w64e", "w64b"), "Huber")

  p <-ggplot(points) +
    geom_segment(aes(x=left, xend=right, y=y, yend=y), lty=points$lty,
                 colour=points$cols, lwd=1)+
    geom_point(aes(x=x, y=y), col="white", shape=20, size=10)+
    geom_point(aes(x=x, y=y), col=points$cols, shape=points$shape, size=2, stroke=1.3)+
    annotate(geom="text", x=c(1, 2, 3), y=rep(2.8, 3),
             label=c("FRET 1", "FRET 2", "FRET 6"), size=3)+
    annotate(geom="text", x=c(1, 2, 3), y=rep(1.8, 3),
             label=c("DESeq2", "WaveQTL", "Huber"), size=3)+
    annotate(geom="text", x=c(4.2, 4.2), y=c(1, 1.5),
             label=c("Naive Bins", "Informed Bins"), size=3)+
    xlim(0.55, 4.5) +
    theme_bw() + theme(panel.grid=element_blank(), axis.title=element_blank(),
                       axis.text=element_blank(),
                       panel.border=element_blank(), axis.ticks=element_blank())
  return(list(p, points))
}

plot_nregion <- function(){

  #FRET
  df = getobj("region_counts.RData")
  names(df)= c("fdr", "total",  "Regions overlapping peaks")
  df$`All regions` = df$total - df$`Regions overlapping peaks`
  ymax = max(df$total)
  dflong = gather(df, "class", "count", -fdr ,-total)
  dflong$class=factor(dflong$class, levels=rev(levels(dflong$class)))
  gridlines = seq(10000, ymax, by=10000)
  h_fret = ggplot(dflong) + geom_area(aes(x=fdr, y = count, fill=class), alpha=0.5) +
    ylab("Number of Regions") + xlab("rFDR Threshold") +
    scale_x_continuous(breaks=seq(0.04, 0.2, by=0.02), limits=c(0.03, 0.2))+
    scale_fill_manual(values=c("darkorange1", "blue"))+
    geom_hline(yintercept = gridlines, lty=3) +
    theme_bw(12) + theme(panel.grid=element_blank(),
                         legend.position=c(0.3, 0.9), legend.title=element_blank())

  #Wellington Bootstrap
  df = getobj("region_counts_well.RData")
  names(df)= c("score", "total",  "Footprints overlapping peaks")
  df=df[df$score >=10,]
  ymax = max(df$total)
  df$`All footprints` = df$total - df$`Footprints overlapping peaks`

  dflong = gather(df, "class", "count", -score ,-total)
  dflong$class=factor(dflong$class, levels=rev(levels(dflong$class)))
  gridlines = seq(10000, ymax, by=10000)
  h_well = ggplot(dflong) + geom_area(aes(x=score, y = count, fill=class), alpha=0.5) +
    ylab("Number of Footprints") + xlab("Score") +
    scale_x_reverse(breaks=seq(10, 100, by=20), limits=c(100, 10))+
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
  names(df)= c("fdr", "total")
  ymax = max(ymax, max(df$total))
  h_wave = ggplot(df) + geom_area(aes(x=fdr, y = total), fill="blue", alpha=0.5) +
    ylab("Number of Peaks") + xlab("FDR Threshold") +
    scale_x_continuous(breaks=c(0.04, 0.06, 0.08, 0.1, 0.15, 0.2), limits=c(0, 0.2))+
    theme_bw(12) + theme(panel.grid=element_blank())


  #Titles and scales
  h_fret = h_fret + ggtitle("(a) FRET")
  h_well = h_well + ggtitle("(b) Wellington-Bootstrap")
  h_huber = h_huber + ggtitle("(b) Huber Fixed-Window Test")
  h_deseq2 = h_deseq2 + ggtitle("(b) DESeq2")
  h_wave = h_wave + ggtitle("(b) WaveQTL")

  ggsave(h_fret, file="~/Dropbox/Thesis/img/fret_results.png", height=5, width=5, units="in", dpi=300)
  ggsave(h_well, file="~/Dropbox/Thesis/img/well_results.png", height=5, width=5, units="in", dpi=300)
  ggsave(h_huber, file="~/Dropbox/Thesis/img/huber_results.png", height=5, width=5, units="in", dpi=300)
  ggsave(h_deseq2, file="~/Dropbox/Thesis/img/deseq2_results.png", height=5, width=5, units="in", dpi=300)
  ggsave(h_wave, file="~/Dropbox/Thesis/img/wave_results.png", height=5, width=5, units="in", dpi=300)
}

#'@export
thin_qqplot <- function(pvals, thin=c(0.25, 100), shade=TRUE, truncate=Inf){
  n <- length(pvals)
  pvals <- sort(pvals)
  v = qunif(p=seq(0, 1, length.out = n+1)[2:(n+1)])
  v = -log10(v)
  q = quantile(v, probs=thin[1])
  q_idx = min(which(v <= q))
  thin_idx = unique(ceiling(seq(n-q_idx, n, length.out=thin[2])))
  thin_idx = c(1:(n-q_idx))
  v=v[thin_idx]
  pvals = -log10(pvals[thin_idx])

  #Shading
  if(shade){
    c975 <- sapply(thin_idx, FUN=function(i){
      qbeta(0.975, i, n - i + 1)
    })
    c025 <- sapply(thin_idx, FUN=function(i){
      qbeta(0.025, i, n - i + 1)
    })
    df.shade = data.frame("x"=c(v, rev(v)), "y"=c(-log10(c025), rev(-log10(c975))))
  }

  if(is.finite(truncate) & any(pvals > truncate)){
    shape <- c(1, 2)[as.numeric(pvals > truncate) + 1]
    pvals <- pmin(truncate, pvals)
  }else{
    shape=1
  }


  df <- data.frame("pval"=pvals, "v"=v)
  h <- ggplot(df)
  if(shade) h <-  h + geom_polygon(data=df.shade, aes(x=x, y=y), alpha=0.3, fill="black")
  h <- h +  geom_point(aes(x=v, y=pval), shape=shape) +
    geom_abline(slope=1, intercept=0) +
    xlab(expression(Expected~~-log[10](italic( p )-value)))+
    ylab(expression(Observed~~-log[10](italic( p )-value)))+
    theme_bw(12) + theme(panel.grid=element_blank())
  return(h)
}
