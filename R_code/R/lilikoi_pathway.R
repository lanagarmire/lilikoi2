#' @name lilikoi_pathway
#' @title
#'
#' @param gene.data gene expression dataset
#' @param cpd.data metabolite compounds dataset
#' @param pathway.id pathway id
#' @return Pathview visualization output
#' @export
#' @examples lilikoi_pathview(gene.data, cpd.data, pathway.id)

library(pathview)

lilikoi_pathview <- function(gene.data=NULL, cpd.data=NULL, pathway.id,
                             species="hsa", out.suffix=NULL, keys.align="y",
                             kegg.native=T, key.pos) {

  pv.out <- pathview(gene.data, cpd.data,
                     pathway.id, species, out.suffix,
                     keys.align, kegg.native, key.pos = demo.paths$kpos1[i])

  return(pv.out)
}

