## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  install.packages("lilikoi")
#  library(lilikoi)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  dt <- lilikoi.Loaddata(file=system.file("extdata", "plasma_breast_cancer.csv", package = "lilikoi"))
#  Metadata <- dt$Metadata
#  dataSet <- dt$dataSet

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  convertResults=lilikoi.MetaTOpathway('name')
#  Metabolite_pathway_table = convertResults$table
#  head(Metabolite_pathway_table)

