
#' A PDSfun Function
#'
#' This function allows you to compute Pathway Desregulation Score deriving
#' make sure that you have the below database for the metabolites and pathway list:
#' meta_path.RData
#' @param qvec This is the Metabolite_pathway_table from MetaTOpathway function. This table includes the metabolites ids and the its corssponding hmdb ids
#' @keywords PDS
#' @export
#' @examples PDSmatrix= PDSfun(Metabolite_pathway_table)
#' PDSfun(Metabolite_pathway_table)
#'
#'

PDSfun<-function(qvec){

  require(pathifier)
  #load("../data/meta_path.RData",.GlobalEnv);
  #load(paste0(getwd(),"/","data/","meta_path.RData"), .GlobalEnv);

  phe=(Metadata$Label) %>% as.numeric  %>% -1

  newData1=qvec %>% filter(pathway!='NA')%>% select(Query,HMDB)
  newData=Metadata[,t(newData1['Query'])]
  colnames(newData)=t(newData1['HMDB'])
  newData=t(newData)

  PDS<-quantify_pathways_deregulation(as.matrix(newData), row.names(newData), metabolites.list,
                                      pathway.list,as.logical(phe), attempts = 5, min_exp=0, min_std=0)

  qpdmat <<- matrix(as.data.frame(PDS$scores), nrow=length(names(PDS$scores)), byrow=TRUE)
  #qpdmat <<- data.frame(lapply(PDS$scores,function (x) {as.numeric(unlist(x))}),check.names=F)
  colnames(qpdmat) <<- colnames(newData)
  rownames(qpdmat) <<- names(PDS$scores)
  mode(qpdmat) <- "numeric"
  return(qpdmat)
}
