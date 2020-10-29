#' lilikoi.KEGGplot
#'
#' Visualizes selected pathways based on their metabolites expression data.
#'
#' @param metamat metabolite expression data matrix
#' @param sampleinfo  is a vector of sample group, with element names as sample IDs.
#' @param grouporder grouporder is a vector with 2 elements,
#' the first element is the reference group name, like 'Normal', the second one is the experimental group name like 'Cancer'.
#' @param pathid character variable, Pathway ID, usually 5 digits.
#' @param specie character, scientific name of the targeted species.
#' @param filesuffix output file suffix
#' @param Metabolite_pathway_table Metabolites mapping table
#' @return Pathview visualization output
#' @import limma
#' @importFrom plyr ddply
#' @importFrom stats model.matrix
#' @export
#' @examples
#' \donttest{
#' dt = lilikoi.Loaddata(file=system.file("extdata","plasma_breast_cancer.csv", package = "lilikoi"))
#' Metadata <- dt$Metadata
#' dataSet <- dt$dataSet
#' convertResults=lilikoi.MetaTOpathway('name')
#' Metabolite_pathway_table = convertResults$table
#'
#' metamat <- Metadata[, -1]
#' sampleinfo <- Metadata$Label
#' names(sampleinfo) <- rownames(Metadata)
#' grouporder <- unique(Metadata$Label)
#' options(bitmapType='cairo')
#' lilikoi.KEGGplot(metamat = metamat, sampleinfo = sampleinfo, grouporder = grouporder,
#' pathid = '00250', specie = 'hsa',filesuffix = 'GSE16873',Metabolite_pathway_table = Metabolite_pathway_table)
#' }

lilikoi.KEGGplot <- function(metamat, sampleinfo, grouporder, pathid = '00250', specie = 'hsa',
                             filesuffix = 'GSE16873',
                             Metabolite_pathway_table = Metabolite_pathway_table){

  meta_table <- Metabolite_pathway_table$table[, c("Query", "KEGG", "pathway")]

  for(i in 1:ncol(meta_table)){

    meta_table[,i] <- as.character(meta_table[,i])

  }

  keggmets <- subset(meta_table, KEGG != '')

  sharedmets <- intersect(colnames(metamat), keggmets$Query)
  keggmets <- subset(keggmets, Query %in% sharedmets)
  metamat <- metamat[,keggmets$Query]
  metamat <- metamat[names(sampleinfo),]
  metamat <- t(metamat)
  metamat <- metamat[keggmets$Query,]


  mergedat <- as.data.frame(metamat)
  mergedat$KEGG <- keggmets$KEGG
  mergedat <- mergedat[c('KEGG', colnames(metamat))]

  mergekegg <- function(sub){

    keggname <- unique(sub[,1])
    subvals <- sub[-1]
    subval <- colMeans(subvals)
    subval <- as.data.frame(subval, stringsAsFactors = FALSE)
    subval <- t(subval)
    subval <- as.data.frame(subval)
    samplenames <- names(subval)
    subval$KEGG <- keggname
    subval <- subval[c('KEGG', samplenames)]
    row.names(subval) <- 1:nrow(subval)
    return(subval)

  }

  mergedat <- ddply(.data = mergedat, .variables = c('KEGG'), .fun = mergekegg)

  row.names(mergedat) <- mergedat$KEGG
  mergedat <- mergedat[-1]
  mergedat <- t(t(mergedat))
  rm(metamat)

  pd <- as.data.frame(sampleinfo, stringsAsFactors = FALSE)
  pd$samplename <- row.names(pd)
  pd$samplegroup <- pd$sampleinfo
  pd <- pd[-1]
  row.names(pd) <- NULL
  pd <- subset(pd, samplegroup %in% grouporder)
  pd$samplegroup <- factor(x = pd$samplegroup, levels = grouporder, ordered = TRUE)

  mergedat <- mergedat[,pd$samplename]

  design <- model.matrix(~ samplegroup, data = pd)

  fit1 <- lmFit(mergedat, design)
  fit2 <- eBayes(fit1)

  logfcres <- topTable(fit2, coef = 2, n = nrow(fit2))
  logfcres <- data.frame(samplename = row.names(logfcres), logFC = logfcres$logFC,
                         stringsAsFactors = FALSE)
  row.names(logfcres) <- logfcres$samplename
  logfcres <- logfcres[-1]
  logfcres <- t(t(logfcres))


  pv.out <- pathview(cpd.data = logfcres, pathway.id = pathid, species = specie, out.suffix = filesuffix)

  print(pv.out)

  return(pv.out)


}

