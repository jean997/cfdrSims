#Coppied from GWAS Tools
#' Get an object from an .RData file
#' @export
getobj <- function (Rdata){
  objname <- load(Rdata)
  if (length(objname) > 1) {
    warning(paste("Multiple objects stored in file", Rdata,
                  "\nReturning only the first object"))
  }
  return(get(objname))
}
