
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
  if(length(ix_pk > 0)) wintest_res[1, 2] = format(peak_stats$HuberQ_05[ix_pk], digits=2)
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
    wintest_res[2, 2] = format(peak_stats$qvalue[ix2], digits=2)
    wintest_res[3, 2] = format(peak_stats$padj[ix2], digits=2)
  }
  rm(peak_stats)

  #Wellington Stats
  wellington = read_delim("Wellington_Footprints/all_footprints.sorted.bed", delim="\t", col_names=FALSE)
  wellington = wellington[wellington$X1==chr,]
  wl = Intervals(wellington[, c(2, 3)])
  ix_wl = unlist(interval_overlap(myI, wl))
  if(length(ix_wl) > 0) wintest_res[4, 1] = paste0(wellington[ix_wl, c(5, 7)], collapse=" ")


  wintest_res = data.frame(wintest_res)
  names(wintest_res)=c("Hotspot", "Peak")
  wintest_res$Test=c("Huber", "DESeq2-Q", "DESeq2 -QIF", "Wellington")
  wintest_res = wintest_res[, c("Test", "Hotspot", "Peak")]

  #Read data
  cmd = paste0("awk '{if($2==", ix_hs, "){print}}' hotspot_dat/hs_", chr, ".txt > temp.dat")
  system(cmd)
  dat = read_delim("temp.dat", delim=" ", col_names=FALSE)
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
  fdr_levels=c(0.5, 0.2, 0.1, 0.08, 0.05, 0.04)
  if(length(ix_fret_chr)==0){
    cat("No FRET statistics in this region. Recalculating\n")
    ss = huber_helper(as.matrix(dat)[, -c(1, 2)], X$X2)
    ss = ss[1,]/(ss[2,] + 0.05)
    stats = ksmooth_0(x=dat$pos, y=ss, bandwidth=50)
    stats.data = data.frame(cbind(dat$pos, stats))[keep,]
    names(stats.data)=c("pos", "stat")
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
    stat_at_fdr = data.frame(cbind( rep(fdr_levels, each=2), c(-1, 1)*stat_at_fdr))
    names(stat_at_fdr) = c("fdr", "stat")
  }
  dat=dat[keep,]
  datlong = gather(dat, "sample", "count", -pos, -win)
  datlong$Sensitve = factor(X$X2[match(datlong$sample, X$nn)])

  dataplot = ggplot(datlong) + geom_line(aes(x=pos, y=count, group=sample, color=Sensitve)) +
        theme_bw() + xlab("Position") + ylab("DNase 1 Sensitivity")

  statplot = ggplot(stat.data) + geom_line(aes(x=pos,  y=stat)) + xlab("Statistic") + ylab("Position") +
  if(length(ix_fret_chr)==1) statplot = statplot + geom_hline(data=stat_at_fdr, aes(yintercept=stat, col=fdr))
  statplot = statplot + geom_hline(yintercept = c(-1, 1)*0.9, col="red") + theme_bw()

  return(list("dataplot"=dataplot, "statplot"=statplot, "wintest_res"=wintest_res))

}
