
#' Plot Simulations
#'@description Plot simulations
#'@param summ.files List of summary files (see collect.R)
#'@param lty.legend Legend labels for line type
#'@param ltys Line type for each summary file
#'@param ymax Optional max height of FDR y axis
#'@param levels Levels in summ.files
#' @return nothing
#'@export
make_plot <- function(tot.rates, names, ltys, cols, shapes,
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

  tpr$type <- fdp$type <- as.factor(names)
  names(tpr) <- names(fdp) <- c(levels, "type")
  tprlong <- gather(tpr, "level", "tpr", -type, convert=TRUE)
  fdplong <- gather(fdp, "level", "fdp", -type, convert=TRUE)

  #tprlong$type <- relevel(tprlong$type, names)
  #fdplong$type <- relevel(fdplong$type, names)

  tprplot <- ggplot(tprlong) + geom_line(aes(x=level, y=tpr, group=type, col=type, lty=type), lwd=2) +
    geom_point(aes(x=level, y=tpr), size=12, colour="white", shape=20) +
    geom_point(aes(x=level, y=tpr, col=type, shape=type), size=4) +
    theme_bw() + labs(x="Target FDR", y="Average True Positive Rate") +
    scale_color_manual(values=cols) +
    scale_linetype_manual(values=ltys) +
    scale_shape_manual(values=shapes) +
    scale_x_continuous(breaks=levels, limits=c(0, max(levels))) +
    scale_y_continuous(limits = c(0, max(tot.rates$avg.tpr[ix,]))) +
    theme(panel.grid=element_blank())

  if(max(tot.rates$avg.fdp[ix,]) < max(levels)) ymx <- max(levels)
    else ymx <- max(tot.rates$avg.fdp[ix,])
  fdpplot <- ggplot(fdplong) + geom_line(aes(x=level, y=fdp, group=type, col=type, lty=type), lwd=2) +
    geom_point(aes(x=level, y=fdp), size=12, colour="white", shape=20) +
    geom_point(aes(x=level, y=fdp, col=type, shape=type), size=4) +
    geom_abline(slope=1, intercept=0) +
    theme_bw() + labs(x="Target FDR", y="Average False Discovery Proportion") +
    scale_color_manual(values=cols) +
    scale_linetype_manual(values=ltys) +
    scale_shape_manual(values=shapes) +
    scale_x_continuous(breaks=levels, limits=c(0, max(levels))) +
    scale_y_continuous(breaks=levels, limits = c(0, ymx)) +
    theme(panel.grid=element_blank())
  return(list("tprplot"=tprplot, "fdpplot"=fdpplot))
}
