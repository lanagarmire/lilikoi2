#' A PDSfun Function
#'
#' This function allows you to compute Pathway Desregulation Score deriving
#' make sure that you have the below database for the metabolites and pathway list:
#' meta_path.RData
#' @param qvec This is the Metabolite_pathway_table from MetaTOpathway function. This table includes the metabolites ids and the its corssponding hmdb ids
#' @keywords PDS
#' @import dplyr pathifier
#' @return A large matrix of the pathway deregulation scores for each pathway in different samples.
#' @export
#' @examples
#' \donttest{
#' dt <- lilikoi.Loaddata(file=system.file("extdata",
#'   "plasma_breast_cancer.csv", package = "lilikoi"))
#' Metadata <- dt$Metadata
#' dataSet <- dt$dataSet
#' convertResults=lilikoi.MetaTOpathway('name')
#' Metabolite_pathway_table = convertResults$table
#' # PDSmatrix= lilikoi.PDSfun(Metabolite_pathway_table)
#' }

lilikoi.PDSfun<-function(qvec){

  Metadata$Label <- as.factor(Metadata$Label)
  phe=(Metadata$Label) %>% as.numeric  %>% -1

  newData1=qvec %>% filter(pathway!='NA')%>% select(Query,HMDB)
  newData=Metadata[,t(newData1['Query'])]
  colnames(newData)=t(newData1['HMDB'])
  newData=t(newData)

  PDS<-quantify_pathways_deregulation(as.matrix(newData), row.names(newData), metabolites.list,
                                      pathway.list,as.logical(phe), attempts = 5, min_exp=0, min_std=0)

  qpdmat <- matrix(as.data.frame(PDS$scores), nrow=length(names(PDS$scores)), byrow=TRUE)
  #qpdmat <<- data.frame(lapply(PDS$scores,function (x) {as.numeric(unlist(x))}),check.names=F)
  colnames(qpdmat) <- colnames(newData)
  rownames(qpdmat) <- names(PDS$scores)
  mode(qpdmat) <- "numeric"
  return(qpdmat)
}
