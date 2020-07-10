#' A Loaddata Function
#'
#' This function allows you to load your metabolomics data.
#' @param filename file name.
#' @keywords Load
#' @importFrom utils read.csv head
#' @return A data frame named Metadata.
#' @export
#' @examples
#' \donttest{
#'  lilikoi.Loaddata(file=system.file("extdata", "plasma_breast_cancer.csv", package = "lilikoi"))
#' }

lilikoi.Loaddata<-function(filename){

  Metadata <- read.csv(file=filename,check.names=F,row.names=1)
  dataSet <- list()
  dataSet$cmpd<-(colnames(Metadata)[-1])
  newList <- list("Metadata"=Metadata, "dataSet"=dataSet)
  return(newList)
}


