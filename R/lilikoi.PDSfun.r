#' A PDSfun Function
#'
#' This function allows you to compute Pathway Desregulation Score deriving
#' make sure that you have the below database for the metabolites and pathway list:
#' meta_path.RData
#' @param qvec This is the Metabolite_pathway_table from MetaTOpathway function. This table includes the metabolites ids and the its corssponding hmdb ids
#' @keywords PDS
#' @import dplyr princurve
#' @return A large matrix of the pathway deregulation scores for each pathway in different samples.
#' @references Nygård, S., Lingjærde, O.C., Caldas, C. et al.
#' PathTracer: High-sensitivity detection of differential pathway activity in tumours.
#' Sci Rep 9, 16332 (2019). https://doi.org/10.1038/s41598-019-52529-3
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
  phe = (Metadata$Label) %>% as.numeric %>% -1
  newData1 = qvec %>% filter(pathway != "NA") %>% select(Query,HMDB)
  newData = Metadata[, t(newData1["Query"])]
  colnames(newData) = t(newData1["HMDB"])
  newData = t(newData)

  ##############
  ## Modified from pathtracer compute.pts()
  ## https://github.com/staaln/pathtracer/blob/master/R/compute.pts.R
  ## Citation:
  ## Nygård, S., Lingjærde, O.C., Caldas, C. et al.
  ## PathTracer: High-sensitivity detection of differential pathway activity in tumours.
  ## Sci Rep 9, 16332 (2019). https://doi.org/10.1038/s41598-019-52529-3
  ##############
  compute.pts = function(dat, reference, ncomp=4) {
    # Check reference argument
    if (mode(reference) == "logical") {
      if (length(reference) != ncol(dat)) {
        stop(paste("Boolean argument has wrong length: reference"))
      } else if (!any(reference)) {
        stop(paste("Reference population must have at least one member"))
      } else if (all(reference)) {
        stop(paste("Non-reference population must have at least one member"))
      }
    } else if (mode(reference) == "numeric") {
      reference = (1:ncol(dat)) %in% reference
    } else {
      stop("Argument must be logical vector or numeric vector: reference")
    }

    # Perform PCA and keep requested number of components
    m = min(ncomp, min(dim(dat)))
    print(dim(dat))
    #print(dat)
    udv = svd(scale(t(dat), scale=F))
    dat.pca = udv$u[,1:m] %*% diag(udv$d[1:m])  # samples x m

    # Compute principal curve and PDS/QDS scores
    res = principal_curve(dat.pca, start=rbind(dat.pca[reference,], dat.pca[!reference,]))
    pds = res$lambda
    pds = pds/max(pds)
    return(pds)
  }

  qpdmat = c()
  for (path in pathway.list) {
    metabolite_temp = unname(metabolites.list[path][[1]])
    dat_temp = newData[rownames(newData)[rownames(newData)%in%metabolite_temp],,drop = FALSE ]
    ncomp = 4
    if(dim(dat_temp)[1]>=3){ #mindim=4:201, mindim=3:252, mindim=2:324
      pds_temp = t(as.matrix(compute.pts(dat_temp,phe,ncomp)))
      rownames(pds_temp) = path
      colnames(pds_temp) = colnames(newData)
      qpdmat = rbind(qpdmat,pds_temp)
    }
  }

  mode(qpdmat) <- "numeric"
  return(qpdmat)
}
