## Lilikoi is a novel tool for personalized pathway analysis of metabolomics data.

Previously we developed Lilikoi, a personalized pathway-based method to classify diseases using metabolomics data. Given the new trends of computation in the metabolomics field, here we report the next version of Lilikoi as a significant upgrade. The new Lilikoi v2 R package has implemented a deep-learning method for classification, in addition to popular machine learning methods. It also has several new modules, including the most significant addition of prognosis prediction, implemented by Cox-PH model and the deep-learning based Cox-nnet model. Additionally, Lilikoi v2 supports data preprocessing, exploratory analysis, pathway visualization and metabolite-pathway regression. In summary, Lilikoi v2 is a modern, comprehensive package to enable metabolomics analysis in R programming environment.

## Installation

```
install.packages("lilikoi")

# Or for the latest dev version:
devtools::install_github("lanagarmire/lilikoi2")
```

## Example

```
# library(lilikoi)

dt <- lilikoi.Loaddata(file=system.file("extdata", "plasma_breast_cancer.csv", package = "lilikoi"))
Metadata <- dt$Metadata
dataSet <- dt$dataSet

# Transform the metabolite names to the HMDB ids using Lilikoi MetaTOpathway function
convertResults=lilikoi.MetaTOpathway('name')
Metabolite_pathway_table = convertResults$table
head(Metabolite_pathway_table)

# Transform metabolites into pathway using pathtracer algorithm
PDSmatrix=lilikoi.PDSfun(Metabolite_pathway_table)

# Select the most signficant pathway related to phenotype.
selected_Pathways_Weka= lilikoi.featuresSelection(PDSmatrix,threshold= 0.50,method="gain")

# Machine learning
lilikoi.machine_learning(MLmatrix = Metadata, measurementLabels = Metadata$Label,
                              significantPathways = 0,
                              trainportion = 0.8, cvnum = 10, dlround=50,nrun=10, Rpart=TRUE,
                              LDA=TRUE,SVM=TRUE,RF=TRUE,GBM=TRUE,PAM=FALSE,LOG=TRUE,DL=TRUE)
                              
# Prognosis model
lilikoi.prognosis(event, time, exprdata, percent=percent, alpha=0, nfold=5, method="quantile",
          cvlambda=cvlambda,python.path=NULL,coxnnet=FALSE,coxnnet_method="gradient")
          
# Metabolites-pathway regression
lilikoi.meta_path(PDSmatrix = PDSmatrix, selected_Pathways_Weka = selected_Pathways_Weka, Metabolite_pathway_table = Metabolite_pathway_table, pathway = "Pyruvate Metabolism")

# KEGG plot
lilikoi.KEGGplot(metamat = metamat, sampleinfo = sampleinfo, grouporder = grouporder,
                 pathid = '00250', specie = 'hsa',
                 filesuffix = 'GSE16873', 
                 Metabolite_pathway_table = Metabolite_pathway_table)
```



# Built By

*   Xinying Fang https://github.com/vivid225
*   Yu Liu 
*   Zhijie Ren 
*   Yuheng Du https://github.com/yhdu36
*   Qianhui Huang
*   Fadhl Alakwaa https://github.com/FADHLyemen
*   Sijia Huang https://github.com/scarlettcanny
*   Lana Garmire https://github.com/lanagarmire

# More Examples

*   https://github.com/lanagarmire/lilikoi2/blob/master/Lilikoi2%20User%20Guide.Rmd
