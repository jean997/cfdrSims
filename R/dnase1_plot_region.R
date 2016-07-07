
#Run in directory that has merged_hotspots.bed
#and hotspot_dat/hs_chr*.txt
dnase1_plot_region <- function(chr, strt, stp,
                               bandwidth=50, buffer=300){

  options("scipen"=3)

  hs = read_delim("merged_hotspots.bed", delim="\t", col_names=FALSE)
  #Figure out which hotspot overlaps
  hs = hs[hs$X1==chr,]
  library(intervals)
  hs = Intervals(hs[, 2:3])
  myI = Intervals(c(strt, stp))
  ix_hs = unlist(interval_overlap(myI, hs))
  stopifnot(length(ix_hs)==1)

  #Interval boundaries
  bounds = data.frame(hs[ix_hs,])
  bounds$type="hotspot"

  wintest_res = matrix(nrow=4, ncol=2)
  ####Hotspot Peak
  #Huber
  #Deseq2 Qvalue
  #DESeq2 padj
  #Wellington

  #Huber Statistics
  hotspot_stats = getobj("hotspot_dat/hs_all_tests_s0-0.05.RData")
  hotspot_stats = hotspot_stats[hotspot_stats$chr==chr,]
  ix2 = which(hotspot_stats$Window==ix_hs)
  wintest_res[1, 1] = format(hotspot_stats$HuberQ_05[ix2], digits=2)
  rm(hotspot_stats)

  peak_stats = getobj("hotspot_dat/peak_all_tests_s0-0.05.RData")
  peak_stats = peak_stats[peak_stats$chr==chr, ]
  pk = Intervals(peak_stats[, c("winstart", "winstop")])
  ix_pk = unlist(interval_overlap(myI, pk))
  if(length(ix_pk > 0)){
    wintest_res[1, 2] = paste0( format(peak_stats$HuberQ_05[ix_pk], digits=2), collapse=";")
    bounds = rbind(bounds, data.frame(pk[ix_pk,], type="peak"))
  }
  rm(peak_stats)

  #DESeq2 Statistics
  hotspot_stats = getobj("deseq2_analysis/hs_deseq2.RData")
  hotspot_stats = hotspot_stats[hotspot_stats$chr==chr,]
  ix2 = which(hotspot_stats$winstart==hs[ix_hs,1])
  wintest_res[2, 1] = format(hotspot_stats$qvalue[ix2], digits=2)
  wintest_res[3, 1] = format(hotspot_stats$padj[ix2], digits=2)
  rm(hotspot_stats)

  peak_stats = getobj("deseq2_analysis/peak_deseq2.RData")
  peak_stats = peak_stats[peak_stats$chr==chr,]
  pk = Intervals(peak_stats[, c("winstart", "winstop")])
  ix_pk = unlist(interval_overlap(myI, pk))
  if(length(ix_pk)> 0){
    wintest_res[2, 2] = paste0(format(peak_stats$qvalue[ix_pk], digits=2), collapse=";")
    wintest_res[3, 2] = paste0(format(peak_stats$padj[ix_pk], digits=2), collapse=";")
  }
  rm(peak_stats)

  #Wellington Stats
  wellington = read_delim("Wellington_Footprints/all_footprints.sorted.bed", delim="\t", col_names=FALSE)
  wellington = wellington[wellington$X1==chr,]
  wl = Intervals(wellington[, c(2, 3)])
  ix_wl = unlist(interval_overlap(myI, wl))
  if(length(ix_wl) > 0) {
    rr = round(wellington[ix_wl, 5], digits=2)
    wintest_res[4, 1] = paste(as.matrix(rr), as.matrix(wellington[ix_wl, 7]), collapse=";", sep="")
    bounds = rbind(bounds, data.frame(wl[ix_wl,], type="wellington"))
  }

  wintest_res = data.frame(wintest_res)
  names(wintest_res)=c("Hotspot", "Peak")
  wintest_res$Test=c("Huber", "DESeq2-Q", "DESeq2 -QIF", "Wellington")
  wintest_res = wintest_res[, c("Test", "Hotspot", "Peak")]

  #Read data
  N = runif(n=1, min=1e5, max=1e9)
  cmd = paste0("awk '{if($2==", ix_hs, "){print}}' hotspot_dat/hs_", chr, ".txt > temp", N, ".dat")
  system(cmd)
  dat = read_delim(paste0("temp", N, ".dat"), delim=" ", col_names=FALSE)
  unlink(paste0("temp", N, ".dat"))
  names(dat)[1:2] = c("pos", "win")
  keep = which(dat$pos >= strt -buffer-bandwidth & dat$pos <= stp + buffer+bandwidth)
  dat = dat[keep,]
  #Smooth for plotting purposes
  if(bandwidth > 0){
    dats = apply(dat[, -c(1, 2)], MARGIN=2, FUN=function(y){
      ksmooth_0(x=dat$pos, y=y, bandwidth=bandwidth)
    })
    dat[, 3:27] <- dats
  }
  keep = which(dat$pos >= strt -buffer & dat$pos <= stp + buffer)
  #Read phenotype
  X <- read_delim("lsd1_pheno2.txt", col_names=FALSE, delim=" ")
  X$nn = paste0("X", 3:27)

  #Get Fret results
  fret_ref  = read_delim(paste0("discopony_output/", chr, "/", chr, "_ref.txt"), delim=" ", col_names=FALSE)
  fret_ref = Intervals(fret_ref[, 2:3])
  ix_fret_chr = unlist(interval_overlap(myI, fret_ref))
  fdr_levels=c(0.1, 0.07, 0.06, 0.05, 0.04)
  if(length(ix_fret_chr)==0){
    cat("No FRET statistics in this region. Recalculating\n")
    ss = huber_helper(as.matrix(dat)[, -c(1, 2)], X$X2)
    ss = ss[1,]/(ss[2,] + 0.05)
    stats = ksmooth_0(x=dat$pos, y=ss, bandwidth=50)
    stat.data = data.frame(cbind(dat$pos, stats))[keep,]
    names(stat.data)=c("pos", "stat")
  }else{
    stopifnot(length(ix_fret_chr)==1)
    R = getobj(paste0("discopony_output/", chr, "/", chr, "_maxes.upd.RData"))
    fret_dat=R[[ix_fret_chr]]
    rm(R)
    ix = which(fret_dat$pos <= stp+buffer & fret_dat$pos >=strt-buffer)
    stat.data <- data.frame(cbind( fret_dat$pos, fret_dat$ys))[ix,]
    names(stat.data) = c("pos", "stat")
    R <- getobj("discopony_thresh_fdr.RData")
    ix_fret = which(names(R$Robs)==fret_dat$file)
    stat_at_fdr = approx(y=R$z[,ix_fret], x=R$fdr, xout = fdr_levels)$y
    #stat_at_fdr = data.frame(cbind(rep(fdr_levels, 2), c(stat_at_fdr, -1*stat_at_fdr)))
    #names(stat_at_fdr) = c("fdr", "stat")
  }
  dat=dat[keep,]
  datlong = gather(dat, "sample", "count", -pos, -win)
  datlong$Sensitve = factor(X$X2[match(datlong$sample, X$nn)])

  bounds$y = rep(1, nrow(bounds))
  bounds$y[bounds$type=="hotspot"] = -1
  bounds$y[bounds$type=="peak"] = -2
  bounds$y[bounds$type=="wellington"] = -3
  bounds$color = c("orange", "skyblue", "violet")[-1*bounds$y]

  dataplot = ggplot(datlong) + geom_line(aes(x=pos, y=count, group=sample, color=Sensitve)) +
        theme_bw(18) + xlab("Position") + ylab("DNase 1 Sensitivity") +
        scale_color_manual(values=c("navyblue", "chartreuse3"))+
    geom_rect(aes(xmin=X1, xmax=X2, ymin=y-0.3, ymax=y+0.3), col="black",
              fill=bounds$color, data=bounds, lwd=0.5, alpha=0.5)+
    theme(legend.position="none", panel.grid=element_blank())

  statplot = ggplot(stat.data) + geom_line(aes(x=pos,  y=stat)) +
    geom_hline(yintercept = c(-1, 1)*0.9, col="red") + theme_bw(18) +
    theme(panel.grid=element_blank()) +
    scale_y_continuous(limits=range(stat.data$stat)) +
    xlab("Position") + ylab("Smoothed Statistic")
  if(length(ix_fret_chr)==1){
    yy = rep(fdr_levels, 2)
    xx = c(stat_at_fdr, -1*stat_at_fdr)
    df = data.frame(stat=xx, fdr=format(yy))
    ax = 0.14*(max(dat$pos)-min(dat$pos))
    statplot = statplot +
      geom_text( data=df, aes(y=stat, label=fdr, x=Inf), hjust=-0.2)+
      geom_hline(yintercept = xx, lty=3) +
      theme(plot.margin=unit(c(1, 5, 2, 1), "lines"))
    gt <- ggplot_gtable(ggplot_build(statplot))
    gt$layout$clip[gt$layout$name=="panel"] <- "off"
    h=grid.arrange(rbind(ggplotGrob(dataplot), gt, size="last"),
                   tableGrob(wintest_res, rows=NULL), nrow=2 , heights=c(8, 2))
  }else{
    h=grid.arrange(rbind(ggplotGrob(dataplot), ggplotGrob(statplot), size="last"),
                   tableGrob(wintest_res, rows=NULL), nrow=2 , heights=c(8, 2))
  }

  return(list("plot"=h, "dataplot"=dataplot, "statplot"=statplot, "wintest_res"=wintest_res))

}
