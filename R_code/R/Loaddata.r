#' A Loaddata Function
#'
#' This function allows you to load your metabolomics data.
#' @param Your file name.
#' @keywords Load
#' @export
#' @examples Loaddata("data_format.csv")
#' Loaddata()
Loaddata<-function(filename){

#install all the Dependencies packages
   #list.of.packages <- c("ggplot2", "caret","dplyr ","pathifier","RWeka","infotheo","pROC","reshape2")
   #new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
   #if(length(new.packages)) install.packages(new.packages)

  Metadata <<- read.csv(file=filename,check.names=F,row.names=1)
  dataSet <<- list()
  dataSet$cmpd<<-(colnames(Metadata)[-1])
  head(Metadata)
}
