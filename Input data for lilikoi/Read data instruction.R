
# For the generation of the supplementary figure, I attached the code in this R file.
# For the generation of Figure 2, 3, 5A and 5B, please refer to the Lilikoi2 User Guide at
# https://github.com/lanagarmire/lilikoi2/blob/master/Lilikoi2%20User%20Guide.Rmd

##### Figure 2 #####
## Correspond to lilikoi function: lilikoi.machine_learning
# Three csv files for each bar graph: performance_training_list.csv, performance_training_opt_list.csv,
# performance_testing_list.csv
# There are 10 simulations in each csv. Within each simulation, the rownames are
# "AUC","Sensitivity","Specificity","Accuracy","Precision","Recall","F1_Statistic","Balanced_Accuracy".
# The final output in the paper are the average of these 10 simulations.

##### Figure 3 #####
## Correspond to lilikoi function: lilikoi.prognosis
tmp= read.table("input data for Fig 3 jcevent.csv", quote="\"", comment.char="")
jcevent = tmp$V1

tmp= read.table("input data for Fig 3 jctime.csv", quote="\"", comment.char="")
jctime = tmp$V1

exprdata_tumor <- read.table("input data for Fig 3 exprdata_tumor.csv", row.names = 1, sep = ",")
exprdata_tumor <- as.matrix(exprdata_tumor)

##### Figure 5A and 5B #####
## Correspond to lilikoi function: lilikoi.meta_path
edgedat <- read.csv("input data for Fig 5A edgedat.csv")
nodedat <- read.csv("input data for Fig 5A edgedat.csv")

Fig5B <- read.csv("input data for Fig 5B.csv")

##### Supplementary figure S1 #####
FigS1 <- read.csv("input data for supplementary Fig S1.csv")
