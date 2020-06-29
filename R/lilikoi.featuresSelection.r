#' A featuresSelection Function
#'
#' This function allows you to reduce the pathway diemsion using xxxx
#' @param PDSmatrix from PDSfun function
#' @param threshold to select the top pathways
#' @param method information gain ("info") or gain ratio ("gain")
#' @import caret RWeka infotheo ggplot2
#' @importFrom stats reorder
#' @keywords features selection
#' @return A list of top metabolites or pathways.
#' @export
#' @examples
#' \donttest{
#' dt <- Loaddata(file=system.file("extdata", "plasma_breast_cancer.csv", package = "lilikoi2"))
#' Metadata <- dt$Metadata
#' dataSet <- dt$dataSet
#' Metabolite_pathway_table=MetaTOpathway('name')
#' PDSmatrix= PDSfun(Metabolite_pathway_table)
#' selected_Pathways_Weka= lilikoi.featuresSelection(PDSmatrix,threshold= 0.54,method="gain")
#' }
#'
#'
lilikoi.featuresSelection <- function(PDSmatrix,threshold= 0.5,method="info"){

  pds_matrix=(as.data.frame(cbind(t(PDSmatrix),Label=Metadata$Label)))
  #head(pds_matrix)
  set.seed(2000)
  training_ID <- createDataPartition(pds_matrix$Label, p = .8,list = FALSE,times = 1)
  training_diagnosis<-pds_matrix[training_ID,]
  #head(training_diagnosis)
  if (method=="info"){
    InfoGainAttributeEval(as.logical(training_diagnosis$Label-1) ~ . , data = training_diagnosis)->infogainfeatures
    selected_pathways<-names(infogainfeatures[infogainfeatures>threshold])}
  else{
    GainRatioAttributeEval(as.logical(training_diagnosis$Label-1) ~ . , data = training_diagnosis)->infogainfeatures
    selected_pathways<-names(infogainfeatures[infogainfeatures>threshold])}


  info.paireddiagnosis.R<-discretize(training_diagnosis[,selected_pathways])
  info.paireddiagnosis.R<-cbind(info.paireddiagnosis.R,as.numeric(as.matrix(training_diagnosis[,ncol(training_diagnosis)])))
  I.R <- mutinformation(info.paireddiagnosis.R,method= "emp")
  I.R.paireddiagnosis<-I.R[,ncol(I.R)]
  #I.R.paireddiagnosis
  theTable <- within(as.data.frame(I.R.paireddiagnosis),
                     I.R.paireddiagnosis <- as.numeric(I.R.paireddiagnosis))
  #theTable
  theTable<-cbind(row.names(theTable),theTable)
  theTable<-theTable[-ncol(I.R),]
  colnames(theTable)[1]<-c("name")
  theTable <- transform(theTable,
                        name = reorder(name,order(I.R.paireddiagnosis, decreasing = TRUE)))
  #theTable

  p <- ggplot(theTable, aes(name, I.R.paireddiagnosis)) + geom_col() + xlab(NULL) +
    ylab(NULL)

  p + theme(axis.text.x = element_text(angle = 90))
  p + coord_flip()
  q <- p + aes(stringr::str_wrap(name, 20), I.R.paireddiagnosis) + ylab("Mutual information") +
    xlab("Pathways")
  plot(q + coord_flip())
  #legend("topright",legend=names(I.R.paireddiagnosis[order(I.R.paireddiagnosis,decreasing=TRUE)])[-1], border=FALSE, cex=0.7)

  return(selected_pathways)
}
