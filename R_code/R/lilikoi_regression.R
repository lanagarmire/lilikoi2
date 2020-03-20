#' Metabolite-pathway regression
#'
#' Perform single variate linear regression between selected pathways and each of their metabolites.
#'
#' @param input Pathway name
#' @param PDSmatrix Pathway deregulation score matrix
#' @param selected_Pathways_Weka Selected top pathways from the featureSelection function
#' @param Metabolite_pathway_table Metabolites mapping table
#' @return A list of regression summaries.
#' @export
#' @examples lilikoi_regression(input="Aminoacyl-tRNA Biosynthesis", PDSmatrix, selected_Pathways_Weka, Metabolite_pathway_table)

lilikoi_regression <- function(input="Aminoacyl-tRNA Biosynthesis", PDSmatrix, selected_Pathways_Weka, Metabolite_pathway_table){
  tPDSmatrix <- t(PDSmatrix)
  PDSmatrix_pathway <- tPDSmatrix[, selected_Pathways_Weka]

  # Metabolites to hmdbid
  metapath_table <- Metabolite_pathway_table[, c(1,3,5)]
  colnm <- as.data.frame(colnames(Metadata)) # colnames of metadata
  colnm <- colnm[2:228,, drop = F]
  colnames(colnm) <- "Query"


  metalist <- lilikoi:::metabolites.list

  l <- metalist[input]

  hmdb <- l[[1]]

  HMDB_inter <- intersect(hmdb, metapath_table$HMDB)
  KEGG_inter <- intersect(hmdb, metapath_table$KEGG)

  filtered_meta_name <- metapath_table[metapath_table$HMDB %in% HMDB_inter | metapath_table$KEGG %in% KEGG_inter,]
  query <- as.character(filtered_meta_name$Query)

  out <- NULL
  for (i in 1:length(query)){
    filtered_meta <- as.data.frame(Metadata[,query[i]], drop=F)
    colnames(filtered_meta) <- query[i]

    filtered_path = as.data.frame(PDSmatrix_pathway[, input], drop=F)

    filtered_meta$res <- filtered_path[,1]

    fit = lm(res ~ ., data = filtered_meta)

    out[[i]] <- fit

  }

  lapply(out, summary)

}


