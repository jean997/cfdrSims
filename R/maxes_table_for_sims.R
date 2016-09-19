#'@export
maxes_table_for_sims <- function(smoothed.stats, pos, zmin, z0=0.3*zmin){

  stopifnot(length(z0)==length(zmin))
  s <- length(zmin)
  stopifnot(s %in% c(1, 2))

  #Find peak heights
  #Intervals defined by z0:
  max1 <- mxlist(smoothed.stats[,1], z0, zmin, return.ix = TRUE)
  if(s==1) max1 <- max1[max1$mx > zmin,]
  else max1 <- max1[max1$mx > zmin[1] | max1$mx < zmin[2],]
  if(nrow(max1) == 0){
    cat("No clusters exceed ", zmin, "\n")
  }
  #Convert index to position
  max1$pos <- pos[max1$ix]
  max1$iv1 <- pos[max1$iv1]
  max1$iv2 <- pos[max1$iv2]
  max1$chr <- "chr1"
  max1 <- max1[, c("mx", "chr", "pos","iv1", "iv2")]
  if(ncol(smoothed.stats)==1){
    R <- list("max1"=max1,
              "z0"=z0, "zmin"=zmin)
    return(R)
  }


  max.perm <- apply(smoothed.stats[,-1], MARGIN=2, FUN=function(yys){
      mxlist(yys, z0, zmin, return.ix = TRUE)
    })
  max.perm <- do.call(rbind, max.perm)

  if(nrow(max.perm) > 0){
    #Convert positions
    max.perm$pos <- pos[max.perm$ix]
    max.perm$iv1 <- pos[max.perm$iv1]
    max.perm$iv2 <- pos[max.perm$iv2]
    max.perm$chr <- "chr1"
    max.perm <- max.perm[order(max.perm$pos, decreasing = FALSE), c("mx", "chr", "pos", "iv1", "iv2")]
  }

  R <- list("max1"=max1, "max.perm"=max.perm,
              "z0"=z0, "zmin"=zmin, "n.perm"=ncol(smoothed.stats)-1)
  return(R)
}
