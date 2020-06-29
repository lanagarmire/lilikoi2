#' A imputation function.
#'
#' This function is used to preprocess data via knn imputation.
#'
#' @param inputdata An expression data frame with samples in the rows, metabolites in the columns
#' @param method The method to be used to process data, including
#' @return A KNN imputed dataset with samples in the rows, metabolites in the columns.
#' @importFrom impute impute.knn
#' @export
#' @examples
#' \donttest{
#' dt <- Loaddata(file=system.file("extdata", "plasma_breast_cancer.csv", package = "lilikoi2"))
#' Metadata <- dt$Metadata
#' dataSet <- dt$dataSet
#' lilikoi.preproc_norm(inputdata=Metadata, method="standard")
#' }


lilikoi.preproc_knn <-function(inputdata=Metadata,
                           method=c("knn")){

  if (!is.element(method,
                  c("knn"))
  ) {
    stop("Invalid process method")
  }

  # vals <- inputdata[2:ncol(inputdata)]

  vals <- inputdata

  # KNN Imputation ####
  if (method == "knn"){
    vals <- t(vals)
    vals <- impute.knn(vals)
    vals <- vals$data
    vals <- as.data.frame(vals)
    output <- t(vals)
    output <- as.data.frame(output)
  }


  # Combine label information back to output ####
  # output$Label<- inputdata$Label


  return(output)

}
