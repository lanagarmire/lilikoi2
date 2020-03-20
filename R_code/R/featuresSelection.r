#' A featuresSelection Function
#'
#' This function allows you to reduce the pathway diemsion using xxxx
#' @param PDSmatrix from PDSfun function
#' @keywords features selection
#' @export
#' @examples selected_Pathways_Weka= featuresSelection(PDSmatrix)
#' featuresSelection(PDSmatrix)
#'
#'
#'
featuresSelection <- function(PDSmatrix,threshold= 0.5,method="info"){
  require(caret)
  require(RWeka)
  require(infotheo)
  require(ggplot2)
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
                     I.R.paireddiagnosis <- as.numeric(I.R.paireddiagnosis,
                                                       levels=names(sort(table(I.R.paireddiagnosis),
                                                                         decreasing=TRUE))))
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
  #legend("topright",legend=names(I.R.paireddiagnosis[order(I.R.paireddiagnosis,decreasing=T)])[-1], border=FALSE, cex=0.7)

  return(selected_pathways)
}
