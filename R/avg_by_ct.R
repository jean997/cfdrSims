
#'@title Get average tpr and fpr rates over replicates. Align rates for fits with the same proportion detected.
#'@param tpr.list A list of vectors giving true positive rates
#'@param fct.list A list of vectors giving counts of false positive
#'@return list with fpr, tpr and s.e
#'@export
avg_by_ct <- function(tpr.list, fct.list){
  B <- length(tpr.list)
  n <- max(unlist(fct.list))
  pts <- c(0, 1:n)
  prop.out <- seq(0, 1, length.out=n+1)
  fct.mat <-tpr.mat <-  matrix(nrow=n+1, ncol=B)
  for(i in 1:B){
    if(all(fct.list[[i]]==0)){
      tpr.mat[,i] <- c(mean(tpr.list[[i]]), rep(NA, n))
    }else{
        apprx.tpr <- approx(x=fct.list[[i]], y=tpr.list[[i]],
                        xout=pts, rule=2, ties=max)
        tpr.mat[,i] <- apprx.tpr$y
    }
  }
  m <- rowMeans(tpr.mat, na.rm=TRUE)
  tot.obs <- rowSums(!is.na(tpr.mat))
  var <- (1/(tot.obs-1))*rowSums((tpr.mat-m)^2, na.rm=TRUE)
  return(list("fct"=pts, "tpr"=m, "s.e"=sqrt(var)))
}
