
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
    theme_bw(15) + labs(x="Target FDR", y="Avg. Proportion of Signals Detected") +
    scale_color_manual(values=cols) +
    scale_linetype_manual(values=ltys) +
    scale_shape_manual(values=shapes) + ggtitle("True Discovery Rate") +
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
    theme_bw(15) + labs(x="Target FDR", y="Avg. False Discovery Proportion") +
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


plot_densities <- function(){

  full_dat <- getobj("~/Desktop/Cluster_FDR/dnase1_results/full_results_signif.RData")

  n <- apply(full_dat[, c(20,  24, 23, 25, 26)], MARGIN=2, FUN=sum)
  top <- rep(c(1, 0), c(sum(n[1:3]), sum(n[4:5])))
  full_dat$mean_med <- full_dat$mean_all/full_dat$median_all
  lab1 <- c("FRET rFDR < 0.05",
            "Huber q-value < 0.05",
            "DESeq2 q-value < 0.05")
  lab2 <- c("WaveQTL p-value < 0.001",
            "Wellington-Bootstrap score > 80")

  #Mean/Median Density Plotts
  df <- data.frame("mean_med" = c(full_dat$mean_med[full_dat$fret_q05_ov],
                                  full_dat$mean_med[full_dat$huber_signif],
                                  full_dat$mean_med[full_dat$deseq2_signif],
                                  full_dat$mean_med[full_dat$waveqtl_p3],
                                  full_dat$mean_med[full_dat$well_g80]),
                   "method" = c(rep(c("fret", "huber", "deseq2", "wave", "well"), n)))
  df$top <- top
  df$method <- factor(df$method, levels = c("fret", "huber", "deseq2", "wave", "well"))
  r <- range(df$mean_med[is.finite(df$mean_med)])

  mean_med <- ggplot(df[top==1,]) + geom_density(aes(x=mean_med, fill=method, linetype=method), alpha=0.3) +
    scale_x_log10(limits=r) + xlab("Ratio of Mean to Median Sensitivity") + ylab("Density")+
    scale_fill_manual(labels=lab1, values=brewer.pal(5, name="Set1")[1:3])  +
    scale_linetype_manual(labels=lab1,values=1:3) +
    theme_bw(12) + theme(panel.grid=element_blank(),
                         legend.title=element_blank(),
                         legend.background=element_blank(),
                         legend.position=c(0.7, 0.85))
  ggsave(plot=mean_med, file="~/Dropbox/Thesis/img/mean_med1.png", height=4.5, width=4.5,
         units="in", dpi=300)

  mean_med <- ggplot(df[top==0,]) + geom_density(aes(x=mean_med, fill=method, linetype=method), alpha=0.3) +
    scale_x_log10(limits=r) + xlab("Ratio of Mean to Median Sensitivity") + ylab("Density")+
    scale_fill_manual(labels=lab2, values=brewer.pal(5, name="Set1")[4:5]) +
    scale_linetype_manual(labels=lab2,values=1:2) +
    theme_bw(12) + theme(panel.grid=element_blank(),
                         legend.title=element_blank(),
                         legend.background=element_blank(),
                         legend.position=c(0.7, 0.85))
  ggsave(plot=mean_med, file="~/Dropbox/Thesis/img/mean_med2.png", height=4.5, width=4.5,
         units="in", dpi=300)

  #Fold change density plots
  full_dat$mean_mean <- full_dat$mean_x0/full_dat$median_x1
  df <- data.frame("mean_mean" = c(full_dat$mean_mean[full_dat$fret_q05_ov],
                                  full_dat$mean_mean[full_dat$huber_signif],
                                  full_dat$mean_mean[full_dat$deseq2_signif],
                                  full_dat$mean_mean[full_dat$waveqtl_p3],
                                  full_dat$mean_mean[full_dat$well_g80]),
                   "method" = c(rep(c("fret", "huber", "deseq2", "wave", "well"), n)))
  df$top <- top
  df$method <- factor(df$method, levels = c("fret", "huber", "deseq2", "wave", "well"))
  r <- range(df$mean_mean[is.finite(df$mean_mean)])

  mean_mean <- ggplot(df[top==1,]) + geom_density(aes(x=mean_mean, fill=method, linetype=method), alpha=0.3) +
    scale_x_log10(limits=r) + xlab("Fold Change Between Trait Groups") + ylab("Density")+
    scale_fill_manual(labels=lab1, values=brewer.pal(5, name="Set1")[1:3])  +
    scale_linetype_manual(labels=lab1,values=1:3) +
    geom_vline(xintercept = 1)+
    theme_bw(12) + theme(panel.grid=element_blank(),
                         legend.title=element_blank(),
                         legend.background=element_blank(),
                         legend.position="none")
  ggsave(plot=mean_mean, file="~/Dropbox/Thesis/img/mean_mean1.png", height=4.5, width=4.5,
         units="in", dpi=300)

  mean_mean <- ggplot(df[top==0,]) + geom_density(aes(x=mean_mean, fill=method, linetype=method), alpha=0.3) +
    scale_x_log10() + xlab("Fold Change Between Trait Groups") + ylab("Density")+
    scale_fill_manual(labels=lab2, values=brewer.pal(5, name="Set1")[4:5]) +
    scale_linetype_manual(labels=lab2,values=1:2) +
    geom_vline(xintercept = 1)+
    theme_bw(12) + theme(panel.grid=element_blank(),
                         legend.title=element_blank(),
                         legend.background=element_blank(),
                         legend.position="none")
  ggsave(plot=mean_mean, file="~/Dropbox/Thesis/img/mean_mean2.png", height=4.5, width=4.5,
         units="in", dpi=300)


  #Min mean density plots
  full_dat$min_mean <- pmin(full_dat$mean_x0, full_dat$mean_x1)
  df <- data.frame("min_mean" = c(full_dat$min_mean[full_dat$fret_q05_ov],
                                   full_dat$min_mean[full_dat$huber_signif],
                                   full_dat$min_mean[full_dat$deseq2_signif],
                                   full_dat$min_mean[full_dat$waveqtl_p3],
                                   full_dat$min_mean[full_dat$well_g80]),
                   "method" = c(rep(c("fret", "huber", "deseq2", "wave", "well"), n)))
  df$top <- top
  df$method <- factor(df$method, levels = c("fret", "huber", "deseq2", "wave", "well"))
  r <- range(df$min_mean[df$min_mean > 0])

  min_mean <- ggplot(df[top==1,]) + geom_density(aes(x=min_mean, fill=method, linetype=method), alpha=0.3) +
    scale_x_log10(limits=r) + xlab("Mean of Lower Sensitivity Group") + ylab("Density")+
    scale_fill_manual(labels=lab1, values=brewer.pal(5, name="Set1")[1:3])  +
    scale_linetype_manual(labels=lab1,values=1:3) +
    theme_bw(12) + theme(panel.grid=element_blank(),
                         legend.title=element_blank(),
                         legend.background=element_blank(),
                         legend.position="none")
  ggsave(plot=min_mean, file="~/Dropbox/Thesis/img/min_mean1.png", height=4.5, width=4.5,
         units="in", dpi=300)

  min_mean <- ggplot(df[top==0,]) + geom_density(aes(x=min_mean, fill=method, linetype=method), alpha=0.3) +
    scale_x_log10(limits=r) + xlab("Mean of Lower Sensitivity Group") + ylab("Density")+
    scale_fill_manual(labels=lab2, values=brewer.pal(5, name="Set1")[4:5]) +
    scale_linetype_manual(labels=lab2,values=1:2) +
    theme_bw(12) + theme(panel.grid=element_blank(),
                         legend.title=element_blank(),
                         legend.background=element_blank(),
                         legend.position="none")
  ggsave(plot=min_mean, file="~/Dropbox/Thesis/img/min_mean2.png", height=4.5, width=4.5,
         units="in", dpi=300)

  #Min mean density plots
  full_dat$max_mean <- pmax(full_dat$mean_x0, full_dat$mean_x1)
  df <- data.frame("max_mean" = c(full_dat$max_mean[full_dat$fret_q05_ov],
                                  full_dat$max_mean[full_dat$huber_signif],
                                  full_dat$max_mean[full_dat$deseq2_signif],
                                  full_dat$max_mean[full_dat$waveqtl_p3],
                                  full_dat$max_mean[full_dat$well_g80]),
                   "method" = c(rep(c("fret", "huber", "deseq2", "wave", "well"), n)))
  df$top <- top
  df$method <- factor(df$method, levels = c("fret", "huber", "deseq2", "wave", "well"))
  r <- range(df$max_mean[df$max_mean > 0])

  max_mean <- ggplot(df[top==1,]) + geom_density(aes(x=max_mean, fill=method, linetype=method), alpha=0.3) +
    scale_x_log10(limits=r) + xlab("Mean of Higher Sensitivity Group") + ylab("Density")+
    scale_fill_manual(labels=lab1, values=brewer.pal(5, name="Set1")[1:3])  +
    scale_linetype_manual(labels=lab1,values=1:3) +
    theme_bw(12) + theme(panel.grid=element_blank(),
                         legend.title=element_blank(),
                         legend.background=element_blank(),
                         legend.position="none")
  ggsave(plot=max_mean, file="~/Dropbox/Thesis/img/max_mean1.png", height=4.5, width=4.5,
         units="in", dpi=300)

  max_mean <- ggplot(df[top==0,]) + geom_density(aes(x=max_mean, fill=method, linetype=method), alpha=0.3) +
    scale_x_log10(limits=r) + xlab("Mean of Lower Sensitivity Group") + ylab("Density")+
    scale_fill_manual(labels=lab2, values=brewer.pal(5, name="Set1")[4:5]) +
    scale_linetype_manual(labels=lab2,values=1:2) +
    theme_bw(12) + theme(panel.grid=element_blank(),
                         legend.title=element_blank(),
                         legend.background=element_blank(),
                         legend.position="none")
  ggsave(plot=max_mean, file="~/Dropbox/Thesis/img/max_mean2.png", height=4.5, width=4.5,
         units="in", dpi=300)







}




plot_median_mean <- function(){
  full_dat <- getobj("~/Desktop/Cluster_FDR/dnase1_results/full_results_signif.RData")

  fret_deseq2 <- ggplot(full_dat[full_dat$fret_signif,]) +
      geom_abline(slope=1, intercept=0, linetype=3) +
        geom_point(aes(x=median_all, y=mean_all, col=deseq2_signif), alpha=0.3)+
        labs(x="Median DNase 1 Sensivity",
              y="Mean DNase 1 Sensitivity",
             title="FRET minimum rFDR < 0.05") +
        scale_color_manual(labels=c("DESeq2 q-value > 0.05", "DESeq2 q-value < 0.05"),
                           values=c("turquoise", "darkorange2")) +
      theme_bw(12) + theme(panel.grid=element_blank(),
                           legend.title=element_blank(),
                           legend.position=c(0.8, 0.8))


  deseq2_fret <- ggplot(full_dat[full_dat$deseq2_signif,]) +
    geom_abline(slope=1, intercept=0, linetype=3) +
    geom_point(aes(x=median_all, y=mean_all, col=fret_signif), alpha=0.3)+
    labs(x="Median DNase 1 Sensivity",
         y="Mean DNase 1 Sensitivity",
         title="DESeq2 q-valule < 0.05") +
    scale_color_manual(labels=c("FRET min rFDR > 0.05", "FRET min rFDR < 0.05"),
                       values=c("darkgrey", "darkorange2")) +
    theme_bw(12) + theme(panel.grid=element_blank(),
                         legend.title=element_blank(),
                         legend.background=element_blank(),
                         legend.position=c(0.3, 0.85))
  ggsave(deseq2_fret, file="~/Dropbox/Thesis/img/deseq2_fret.png",
         height=4, width=4, units="in", dpi=300)


  fret_huber <- ggplot(full_dat[full_dat$fret_signif,]) +
    geom_abline(slope=1, intercept=0, linetype=3) +
    geom_point(aes(x=median_all, y=mean_all, col=huber_signif), alpha=0.3)+
    labs(x="Median DNase 1 Sensivity",
         y="Mean DNase 1 Sensitivity",
         title="FRET minimum rFDR < 0.05") +
    scale_color_manual(labels=c("Huber q-value > 0.05", "Huber q-value < 0.05"),
                       values=c("turquoise", "darkorange2")) +
    theme_bw(12) + theme(panel.grid=element_blank(),
                         legend.title=element_blank(),
                         legend.position=c(0.8, 0.8))


  huber_fret <- ggplot(full_dat[full_dat$huber_signif,]) +
    geom_abline(slope=1, intercept=0, linetype=3) +
    geom_point(aes(x=median_all, y=mean_all, col=fret_signif), alpha=0.3)+
    labs(x="Median DNase 1 Sensivity",
         y="Mean DNase 1 Sensitivity",
         title="Huber fixed-window q-valule < 0.05") +
    scale_color_manual(labels=c("FRET min rFDR > 0.05", "FRET min rFDR < 0.05"),
                       values=c("darkgrey", "darkorange2")) +
    theme_bw(12) + theme(panel.grid=element_blank(),
                         legend.title=element_blank(),
                         legend.background=element_blank(),
                         legend.position=c(0.3, 0.85))
  ggsave(huber_fret, file="~/Dropbox/Thesis/img/huber_fret.png",
         height=4, width=4, units="in", dpi=300)


  full_dat$med_mean <- full_dat$median_all/full_dat$mean_all


}

