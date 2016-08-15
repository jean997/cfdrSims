#Functions for calculating false positive rate and true positive rate

tpr.func <- function(x, labels){
    return(sum(x==1 & labels==1)/sum(labels==1))
}
fpr.func <- function(x, labels){
    return(sum(x==1 & labels==0)/sum(labels==0))
}

get_tpr_fpr <- function(sep, labels){
#sep should be p x B labels should be length p
    tpr <- apply(sep, MARGIN=2, FUN=tpr.func, labels=labels)
    fpr <- apply(sep, MARGIN=2, FUN=fpr.func, labels=labels)
		return(list("tpr"=tpr, "fpr"=fpr))
}

#For ROC curves
#'@title Get average tpr and fpr rates over replicates
#'@param tpr.list A list of vectors giving true positive rates
#'@param fpr.list A list of vectors giving false positive rates
#'@param direction Average over fixed fpr (vertical) or fixed tpr (horizontal)
#'@return list with fpr, tpr and s.e
#'@export
avg_by_interp <- function(tpr.list, fpr.list, direction="vertical", npoints=200){
	B <- length(tpr.list)
	if(direction=="vertical"){
		fpr.out <- seq(0, 1, length.out=npoints)
		tpr.mat <- matrix(0, npoints, B)
		for(i in 1:B){
			apprx.tf <- approx(x=fpr.list[[i]], y=tpr.list[[i]], xout=fpr.out, ties=max, rule=2)
			tpr.mat[,i] <- apprx.tf$y
		}
		m <- rowMeans(tpr.mat, na.rm=TRUE)
		tot.obs <- rowSums(!is.na(tpr.mat))
		var <- (1/(tot.obs-1))*rowSums((tpr.mat-m)^2, na.rm=TRUE)
		return(list("fpr"=fpr.out, "tpr"=m, "s.e"=sqrt(var)))
	}else if(direction=="horizontal"){
		tpr.out <- seq(0, 1, length.out=npoints)
		fpr.mat <- matrix(0, npoints, B)
		for(i in 1:B){
			apprx.tf <- approx(x=tpr.list[[i]], y=fpr.list[[i]], xout=tpr.out, ties=max)
			fpr.mat[,i] <- apprx.tf$y
		}
		m <- rowMeans(fpr.mat, na.rm=TRUE)
		tot.obs <- rowSums(!is.na(fpr.mat))
		var <- (1/(tot.obs-1))*rowSums((fpr.mat-m)^2, na.rm=TRUE)
		return(list("fpr"=m, "tpr"=tpr.out, "s.e"=sqrt(var)))
	}
}


