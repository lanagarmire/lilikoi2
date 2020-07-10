#' A Normalization function.
#'
#' This function is used to preprocess data via normalization.
#' It provides three normalization methods: standard normalization, quantile
#' normalization and median fold normalization. The median fold normalization is adapted
#' from http://www.metabolomics-forum.com/index.php?topic=281.0.
#'
#' @param inputdata An expression data frame with samples in the rows, metabolites in the columns
#' @param method The method to be used to process data, including standard normalization (standard), quantile
#' normalization (quantile) and median fold normalization (median).
#' @return A normalized dataset with samples in the rows, metabolites in the columns.
#' @export
#' @importFrom preprocessCore normalize.quantiles
#' @examples
#' \donttest{
#' dt <- lilikoi.Loaddata(file=system.file("extdata",
#'  "plasma_breast_cancer.csv", package = "lilikoi"))
#' Metadata <- dt$Metadata
#' dataSet <- dt$dataSet
#' lilikoi.preproc_norm(inputdata=Metadata, method="standard")
#' }


lilikoi.preproc_norm <-function(inputdata=Metadata,
                           method=c("standard","quantile","median")){

  if (!is.element(method,
                  c("standard","quantile","median","knn"))
  ) {
    stop("Invalid process method")
  }

  vals <- inputdata[2:ncol(inputdata)]

  # vals <- inputdata # Metadata = inputdata

  # Standard normalization ####
  if (method == "standard"){
    standardnorm <- function(colval){

      avg <- mean(colval)
      colsd <- sum((colval - avg)^2)/(length(colval)-1)
      colstandard <- (colval - avg)/sqrt(colsd)

      return(colstandard)
    }
    normvals <- apply(X = vals, MARGIN = 2, FUN = standardnorm)
    output <- as.data.frame(normvals)
  }

  # Quantile normalization ####
  if (method == "quantile"){

    samplenames <- rownames(vals)
    metabolitenames <- colnames(vals)
    output <- normalize.quantiles(as.matrix(t(vals)),copy=F)
    output <- t(output)
    rownames(output) <- samplenames
    colnames(output) <- metabolitenames
  }


  # Median fold normalization ####
  if (method == "median"){
    medianfoldnorm <- function(mat) {
      # Perform median fold change normalisation
      #          X - data set [Variables & Samples]
      medSam <- apply(mat, 1, median)
      medSam[which(medSam==0)] <- 0.0001
      mat <- apply(mat, 2, function(mat, medSam){
        medFDiSmpl <- mat/medSam
        vec<-mat/median(medFDiSmpl)
        return(vec)
      }, medSam)
      return (mat)
    }
    output <- medianfoldnorm(vals)
  }


  # Combine label information back to output ####
  output <- as.data.frame(output)
  output$Label<- as.data.frame(inputdata$Label)


  return(output)

}


