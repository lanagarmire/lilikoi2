#' Metabolite-pathway regression
#'
#' Performs single variate linear regression between selected pathways and each of their metabolites.
#' Output the network plot between pathways and metabolties.
#' @param PDSmatrix Pathway deregulation score matrix
#' @param selected_Pathways_Weka Selected top pathways from the featureSelection function
#' @param Metabolite_pathway_table Metabolites mapping table
#' @importFrom stats lm
#' @importFrom utils type.convert
#' @import scales RCy3
#' @return A bipartite graph of the relationships between pathways and their corresponding metabolites.
#' @export


meta_path <- function(PDSmatrix, selected_Pathways_Weka, Metabolite_pathway_table){

  regression <- function(input, PDSmatrix, selected_Pathways_Weka, Metabolite_pathway_table){
    tPDSmatrix <- t(PDSmatrix)
    PDSmatrix_pathway <- tPDSmatrix[, selected_Pathways_Weka]

    # Metabolites to hmdbid
    metapath_table <- Metabolite_pathway_table[, c(1,3,5)]
    colnm <- as.data.frame(colnames(Metadata)) # colnames of metadata
    colnm <- colnm[2:228,, drop = F]
    colnames(colnm) <- "Query"

    metalist <- metabolites.list

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


  output <- function(reg){
    resbar <- NULL
    nms <- NULL
    for (i in 1:length(reg)){
      i=1
      mat <- as.data.frame(reg[[i]]$coefficients)
      temp1 <- rownames(mat)[2]
      nms <- rbind(nms, temp1)
      temp2 <- reg[[i]]$coefficients[2,1]
      resbar <- rbind(resbar, temp2)
    }

    res <- cbind(nms, resbar)
    temp <- data.frame(res, stringsAsFactors = FALSE)
    temp[] <- lapply(temp, type.convert)

    temp <- temp[order(temp$X2),]

    return(temp)
  }

  # Output nodedat and edgedat
  y <- NULL
  for (i in 1:length(selected_Pathways_Weka)){
    # i = 1
    reg <- regression(input=selected_Pathways_Weka[i], PDSmatrix, selected_Pathways_Weka, Metabolite_pathway_table)
    res <- output(reg)
    path <- rep(selected_Pathways_Weka[i], length(reg))
    res$path <- as.factor(path)

    y <- rbind(y, res)


  }

  edgedat <- y[c(1,3,2)]
  names(edgedat) <- c('source', 'target','weight')
  rownames(edgedat) <- NULL
  edgedat$edge.col <- ifelse(edgedat$weight<=0, "Red", "Blue")

  source <- as.data.frame(edgedat$source,stringsAsFactors=FALSE)
  names(source) <- "nodes"
  target <- as.data.frame(edgedat$target,stringsAsFactors=FALSE)
  names(target) <- "nodes"

  nodedat <- rbind(source, target)
  nodedat <- nodedat %>%
    distinct(nodes)

  ## individual pathway plot
  indres <- subset(y, path==pathway)
  p <- ggplot(data=indres, aes(x=reorder(X1, X2, sum), y=X2)) + geom_bar(stat="identity")
  bipartite.plot = p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10),
                             plot.title = element_text(color="black", hjust = 0.5)) +
    labs(title=pathway, x=NULL, y=NULL)


  ## Bipartite plot
  cytoimage <- function(nodedat, edgedat){


    cytoscapePing()

    cytoscapeVersionInfo ()

    if('weight' %in% names(edgedat)){
      edgedat$absweight <- abs(edgedat$weight)
    }


    names(nodedat)[1] <- 'id'
    names(edgedat)[c(1, 2)] <- c('source', 'target')
    edgedat$interaction <- 1


    #New code#########################################################################################
    improvenodename <- function(nodes = nodedat, edges = edgedat){

      nodeids <- nodes$id

      ids <- gsub(pattern = "\`", replacement = '', x = nodeids, fixed = FALSE)
      ids <- gsub(pattern = "\"", replacement = '', x = ids, fixed = FALSE)
      ids <- gsub(pattern = "\'", replacement = '', x = ids, fixed = FALSE)

      nodes$id <- ids

      edgesources <- edges$source
      edgetargets <- edges$target

      sources <- gsub(pattern = "\`", replacement = '', x = edgesources, fixed = FALSE)
      sources <- gsub(pattern = "\"", replacement = '', x = sources, fixed = FALSE)
      sources <- gsub(pattern = "\'", replacement = '', x = sources, fixed = FALSE)

      edges$source <- sources

      targets <- gsub(pattern = "\`", replacement = '', x = edgetargets, fixed = FALSE)
      targets <- gsub(pattern = "\"", replacement = '', x = targets, fixed = FALSE)
      targets <- gsub(pattern = "\'", replacement = '', x = targets, fixed = FALSE)

      edges$target <- targets

      return(list(nodedat = nodes, edgedat = edges))

    }


    res <- improvenodename(nodes = nodedat, edges = edgedat)

    nodedat <- res$nodedat
    edgedat <- res$edgedat

    rm(res)
    ####################################################################################################




    createNetworkFromDataFrames(nodedat, edgedat, title = 'lilikoinet', collection = 'lilikoi2')


    #Create own style

    stylename <- 'lilikoistyle'


    defaults <- list(NODE_SHAPE = 'ELLIPSE',
                     EDGE_TRANSPARENCY = 255,
                     EDGE_BEND = "0.728545744495502,-0.684997151948455,0.6456513365424503",
                     EDGE_CURVED = TRUE)

    nodelabel <- mapVisualProperty('node label','id','p')

    nodelabelfontsize <- mapVisualProperty('node label font size', 'info', 'd',
                                           c('metabolite', 'pathway'), c(10, 10))

    nodefill <- mapVisualProperty('node fill color', 'info', 'd',
                                  c('metabolite', 'pathway'), c('cyan', 'yellow'))

    nodeheight <- mapVisualProperty('node height', 'info', 'd',
                                    c('metabolite', 'pathway'), c(30, 50))

    nodewidth <- mapVisualProperty('node width', 'info', 'd',
                                   c('metabolite', 'pathway'), c(40, 60))



    createVisualStyle(stylename, defaults, list(nodelabel, nodelabelfontsize, nodefill, nodeheight, nodewidth))


    #Implement the network

    setVisualStyle(stylename)

    lockNodeDimensions(FALSE, style.name = stylename)

    if('edge.col' %in% names(edgedat)){
      unicols <- unique(edgedat$edge.col)

      if(length(unicols) == 2){

        setEdgeColorMapping('edge.col', c('Red', 'Blue'), c('#FF0000', '#0000FF'), 'd', style.name = stylename)

      }else if(length(unicols) > 2){

        setEdgeColorMapping('edge.col', unicols, hue_pal()(length(unicols), 'd', style.name = stylename))

      }

    }

    if('absweight' %in% names(edgedat)){
      setEdgeLineWidthMapping('absweight', range(edgedat$absweight), c(2, 5), 'c', style.name = stylename)
    }


    layoutNetwork('cose')

    exportImage(type = 'PDF')

  }

  nodegraph <- cytoimage(nodedat, edgedat)
  # return(nodegraph)

  returnList <- list()
  returnList$bipartite.plot <- bipartite.plot
  returnList$node.graph <- nodegraph
  return(returnList)
}





