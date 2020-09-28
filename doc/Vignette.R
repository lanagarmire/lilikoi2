## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  install.packages("lilikoi")
#  library(lilikoi)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  dt <- lilikoi.Loaddata(file=system.file("extdata", "plasma_breast_cancer.csv", package = "lilikoi2"))
#  Metadata <- dt$Metadata
#  dataSet <- dt$dataSet

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  convertResults=lilikoi.MetaTOpathway('name')
#  Metabolite_pathway_table = convertResults$table
#  head(Metabolite_pathway_table)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  selected_Pathways_Weka= lilikoi.featuresSelection(PDSmatrix,threshold= 0.54,method="gain")
#  selected_Pathways_Weka

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  # Standard Normalization
#  lilikoi.preproc_norm(inputdata=Metadata, method="standard")
#  lilikoi.preproc_norm(inputdata=Metadata, method="quantile")
#  lilikoi.preproc_norm(inputdata=Metadata, method="median")

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  # KNN Imputation
#  lilikoi.preproc_knn(inputdata=Metadata,method=c("knn"))

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  lilikoi.explr(data, demo.data, pca=TRUE, tsne=FALSE)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  lilikoi.machine_learning(MLmatrix = Metadata, measurementLabels = Metadata$Label,
#                                significantPathways = 0,
#                                trainportion = 0.8, cvnum = 10)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  # Set up prognosis function arguments
#  # Before running Cox-nnet, users need to provide the directory for python3 and the inst file in lilikoi
#  path = path.package('lilikoi', quiet = FALSE) # path = "lilikoi/inst/", use R to run
#  path = file.path(path, 'inst')
#  
#  python.path = "/Library/Frameworks/Python.framework/Versions/3.8/bin/python3"
#  
#  
#  event = jcevent
#  time = jctime
#  percent = NULL
#  exprdata = exprdata_tumor
#  alpha = 0
#  nfold = 5
#  method = "quantile"
#  cvlambda = NULL
#  coxnnet = TRUE
#  coxnnet_method = "gradient"
#  
#  library(reticulate)
#  
#  lilikoi.prognosis(event, time, exprdata, percent=percent, alpha=0, nfold=5, method="quantile",
#            cvlambda=cvlambda,python.path=python.path,path=path,coxnnet=TRUE,coxnnet_method="gradient")
#  

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  metamat <- t(t(Metadata[, -1]))
#  metamat <- log2(metamat)
#  sampleinfo <- Metadata$Label
#  names(sampleinfo) <- rownames(Metadata)
#  grouporder <- unique(Metadata$Label)
#  
#  lilikoi.KEGGplot(metamat = metamat, sampleinfo = sampleinfo, grouporder = grouporder,
#                   pathid = '00250', specie = 'hsa',
#                   filesuffix = 'GSE16873',
#                   Metabolite_pathway_table = Metabolite_pathway_table)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  lilikoi.meta_path(PDSmatrix = PDSmatrix, selected_Pathways_Weka = selected_Pathways_Weka, Metabolite_pathway_table = Metabolite_pathway_table)

