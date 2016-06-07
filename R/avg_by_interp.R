#Functions for calculating false positive rate and true positive rate

tpr.func <- function(x, labels){
    return(sum(x==1 & labels==1)/sum(labels==1))
}
fpr.func <- function(x, labels){
    return(sum(x==1 & labels==0)/sum(labels==0))
}

#For ROC curves
#'@title Get average tpr and fpr rates over replicates
#'@param tpr.list A list of vectors giving true positive rates
#'@param fpr.list A list of vectors giving false positive rates
#'@param npoints Number of points to evaluate
#'@return list with fpr, tpr and s.e
#'@export
avg_by_interp <- function(tpr.list, fpr.list,  npoints=200){
	B <- length(tpr.list)

	fpr.out <- seq(0, 1, length.out=npoints)
	tpr.mat <- matrix(0, npoints, B)
	for(i in 1:B){
		apprx.tf <- approx(x=fpr.list[[i]], y=tpr.list[[i]], xout=fpr.out, ties=max)
		tpr.mat[,i] <- apprx.tf$y
	}
	m <- rowMeans(tpr.mat, na.rm=TRUE)
	tot.obs <- rowSums(!is.na(tpr.mat))
	var <- (1/(tot.obs-1))*rowSums((tpr.mat-m)^2, na.rm=TRUE)
	return(list("fpr"=fpr.out, "tpr"=m, "s.e"=sqrt(var)))
}


