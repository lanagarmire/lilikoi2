#' A preprocessing function.
#' 
#' This function is used to preprocess data via normalization and imputation.
#' It provides three normalization methods: standard normalization, quantile 
#' normalization and median fold normalization. It also enables KNN imputation.
#' 
#' @param inputdata An expression data frame with samples in the rows, metabolites in the columns
#' @param method The method to be used to process data, including   
#' @return a processed dataset with samples in the rows, metabolites in the columns.
#' @export
#' @examples preproc(inputdata=Metadata, method="standard")


lilikoi_preproc <-function(inputdata=Metadata, 
                           method=c("standard","quantile","median","knn")){
  
  if (!is.element(method, 
                  c("standard","quantile","median","knn"))
  ) {
    stop("Invalid process method")
  }
  
  # vals <- inputdata[2:ncol(inputdata)] 
  
  vals <- inputdata # Metadata = inputdata
  
  # Standard normalization ####
  if (method == "standard"){
    standardnorm <- function(colval){
      
      avg <- mean(colval)
      colsd <- sum((colval - avg)^2)/(length(colval)-1)
      colstandard <- (colval - avg)/sqrt(colsd)
      
      return(colstandard)
    }
    vals <- t(vals)
    normvals <- apply(X = vals, MARGIN = 2, FUN = standardnorm)
    normvals <- t(normvals)
    output <- as.data.frame(normvals)
  }
  
  # Quantile normalization ####
  if (method == "quantile"){
    library(preprocessCore)
    
    samplenames <- rownames(vals)
    metabolitenames <- colnames(vals)
    output <- normalize.quantiles((t(vals))) ### It's working on cols or rows?????
    output <- t(output)
    rownames(output) <- samplenames
    colnames(output) <- metabolitenames
  }
  
  
  # Median fold normalization ####
  if (method == "median"){
    vals <- t(vals)
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
    output <- t(output)
  }
  
  
  # KNN Imputation ####
  if (method == "knn"){
    library(impute)
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


