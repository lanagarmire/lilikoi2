##### Figure 2 #####
# Three csv files for each bar graph: performance_training_list.csv, performance_training_opt_list.csv, 
# performance_testing_list.csv
# There are 10 simulations in each csv. Within each simulation, the rownames are 
# "AUC","Sensitivity","Specificity","Accuracy","Precision","Recall","F1_Statistic","Balanced_Accuracy"

##### Figure 3 #####
tmp= read.table("~/Documents/GitHub/lilikoi2/jcevent.csv", quote="\"", comment.char="")
jcevent = tmp$V1

tmp= read.table("~/Documents/GitHub/lilikoi2/jctime.csv", quote="\"", comment.char="")
jctime = tmp$V1

exprdata_tumor <- read.table("exprdata_tumor.csv", row.names = 1, sep = ",")
exprdata_tumor <- as.matrix(exprdata_tumor)

##### Figure 5A and 5B #####
PDSmatrix = read.table("PDSmatrix.csv", row.names = 1, header = TRUE, sep = ",")
tmp = read.csv("selected_Pathways_Weka.csv")
selected_Pathways_Weka = as.character(tmp$x)
Metabolite_pathway_table <- read.csv("Metabolite_pathway_table.csv")

##### Supplementary figure S1 #####
# The function is not within the lilikoi package, but the generation code is the same as the code to generate figure 5B.
library(lilikoi)

dt <- lilikoi.Loaddata(file=system.file("extdata", "plasma_breast_cancer.csv", package = "lilikoi"))
Metadata <- dt$Metadata
dataSet <- dt$dataSet

convertResults=lilikoi.MetaTOpathway('name')
Metabolite_pathway_table = convertResults$table
PDSmatrix= lilikoi.PDSfun(Metabolite_pathway_table)
selected_Pathways_Weka= lilikoi.featuresSelection(PDSmatrix,threshold= 0.54,method="gain")

regression <- function(input, PDSmatrix, selected_Pathways_Weka, Metabolite_pathway_table){
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

# t<- meta_path(input="Alanine, Aspartate And Glutamate Metabolism", PDSmatrix, selected_Pathways_Weka, Metabolite_pathway_table)


## Regression #####
reg1 <- regression(input="Alanine, Aspartate And Glutamate Metabolism",
                   PDSmatrix, selected_Pathways_Weka, Metabolite_pathway_table)
reg2 <- regression(input="Aminoacyl-tRNA Biosynthesis",
                   PDSmatrix, selected_Pathways_Weka, Metabolite_pathway_table)
reg3 <- regression(input="Aspartate Metabolism",
                   PDSmatrix, selected_Pathways_Weka, Metabolite_pathway_table)
reg4 <- regression(input="Biosynthesis Of Amino Acids",
                   PDSmatrix, selected_Pathways_Weka, Metabolite_pathway_table)
reg5 <- regression(input="Canavan Disease",
                   PDSmatrix, selected_Pathways_Weka, Metabolite_pathway_table)
reg6 <- regression(input="Glycine, Serine And Threonine Metabolism",
                   PDSmatrix, selected_Pathways_Weka, Metabolite_pathway_table)
reg7 <- regression(input="Hypoacetylaspartia",
                   PDSmatrix, selected_Pathways_Weka, Metabolite_pathway_table)
reg8 <- regression(input="Metabolic Pathways",
                   PDSmatrix, selected_Pathways_Weka, Metabolite_pathway_table)
reg9 <- regression(input="Nicotinate And Nicotinamide Metabolism",
                   PDSmatrix, selected_Pathways_Weka, Metabolite_pathway_table)
reg10 <- regression(input="Protein Digestion And Absorption",
                    PDSmatrix, selected_Pathways_Weka, Metabolite_pathway_table)


output <- function(reg){
  resbar <- NULL
  nms <- NULL
  for (i in 1:length(reg)){
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

y <- NULL
for (i in 1:length(selected_Pathways_Weka)){
  # i = 1
  reg <- lilikoi_regression(input=selected_Pathways_Weka[i], PDSmatrix, selected_Pathways_Weka, Metabolite_pathway_table)
  res <- output(reg)
  path <- rep(selected_Pathways_Weka[i], length(reg))
  res$path <- as.factor(path)
  
  y <- rbind(y, res)
  
  
}

edgedat <- y[c(1,3,2)]
names(edgedat) <- c('source', 'target','weight')
rownames(edgedat) <- NULL

source <- as.data.frame(edgedat$source,stringsAsFactors=FALSE)
names(source) <- "nodes"
target <- as.data.frame(edgedat$target,stringsAsFactors=FALSE)
names(target) <- "nodes"

library(dplyr)
nodedat <- rbind(source, target)
nodedat <- nodedat %>%
  distinct(nodes)


## Regression plot - bar plot #####
library(ggplot2)
# Basic barplot
p <- ggplot(data=y, aes(x=reorder(X1, X2, sum), y=X2)) + geom_bar(stat="identity")
p + facet_wrap(~path, scales = "free", nrow = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5)) +
  xlab("Metabolites") + ylab(NULL)



