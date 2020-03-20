#' Model Adjustemnt function using clinical factors 
#'
#' This function for adjusted the best performed model using the clinical factors inserted by the user
#' It plots ROC for three models: model1 build using only seelcted pathways, model2 build using clinical factors 
#' and model3 build using selected pathways and the clinical factors  
#' @param PDSmatrix the PDS matrix geneerated using PDSfun function 
#' @param selected_Pathways_Weka selected pathway using WEKA algorithm geneerated from featuresSelection function
#' @param clinical_factors_data the metadata for the samples
#' @param factors which the users want to add to the model 
#' @keywords adjustment
#' @export
#' @examples model_adjustment(result,PDSmatrix,selected_Pathways_Weka,clinical_factors_data,factors=c('Age','Race'))
#' model_adjustment(result,PDSmatrix,selected_Pathways_Weka,clinical_factors_data,factors=c('Age','Race'))
#'
#'
#'

model_adjustment <-function(result,PDSmatrix,selected_Pathways_Weka,clinical_factors_data,factors){
 require(pROC)   
prostate_df=data.frame(t((PDSmatrix[selected_Pathways_Weka,])),Label=Metadata$Label, check.names=T)
colnames(prostate_df)[which(names(prostate_df) == "Label")]='subtype'
irisTrain <- prostate_df[ result$train_inx,]
irisTest  <- prostate_df[-result$train_inx,]
    
best_model=result$models[which.max(result$performance[1,])] # the best model has the high AUC
method=(unlist(best_model)[[1]])
    
    
cartClasses <- predict( best_model, newdata = irisTest,type="prob")
cartClasses1 <- predict( best_model, newdata = irisTest)
cartConfusion=confusionMatrix(data = unlist(cartClasses1), irisTest$subtype)
#ROC_pathway <- roc(predictor=cartClasses[[1]]$Normal,response=irisTest$subtype,levels=rev(levels(irisTest$subtype)))
ROC_pathway <- roc(predictor=as.numeric(unlist(cartClasses[[1]][1])),
                   response=irisTest$subtype,levels=rev(levels(irisTest$subtype)))
#plot(smooth(ROC,method="fitdistr"),print.auc=TRUE,col="green")
smooth_method="binormal" 
 
 #pdf("factors.pdf",width=10,height=10)   
plot(pROC::smooth(ROC_pathway,method=smooth_method),col="black",cex.lab=1.5)
#plot(ROC_pathway,col="black")
par(new=TRUE)
train_index=result$train_inx
    
factor_data=cbind(Label=clinical_factors_data[,2],clinical_factors_data[,factors])
colnames(factor_data)[which(names(factor_data) == "Label")]='subtype'

ROC_factor=createthemodel( factor_data,train_index,method)
    
plot(pROC::smooth(ROC_factor$ROC,method=smooth_method),col="red",cex.lab=1.5)
#plot(ROC_factor,col="red",print.auc=T)
par(new=TRUE)

pathway_factors_data=cbind(factor_data[-1],prostate_df)
colnames(pathway_factors_data)[which(names(pathway_factors_data) == "Label")]='subtype'
    #print(head(pathway_factors_data))
ROC_pathway_factors=createthemodel(pathway_factors_data,train_index,method)
plot(ROC_pathway_factors$ROC,col="blue",cex.lab=1.5)

#plot(ROC_pathway_factors,col="blue") 
    
 legend(0.5, 0.4, legend=c('Selected Pathways','Clinical factors','Selected pathways + clinical factors'), 
 col=c("black", "red","blue"), lty=1:2, cex=1.2)  
 
    #dev.off()
    
    plot(plot(varImp(ROC_factor$model, scale = FALSE,top=20),main=method))
    plot(plot(varImp(ROC_pathway_factors$model, scale = FALSE,top=20),main=method))
# plot the correlation between the pathways and the clinical factors 
    
require(corrplot)
require("Hmisc")
df=select(pathway_factors_data,-subtype)
df=sapply(df,as.numeric)
all_corr=rcorr(as.matrix(df),type="pearson")
    #pdf("corr.pdf",width=20,height=15)
    corrplot(all_corr$r, method='circle',type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45,tl.cex=1,p.mat = all_corr$P,sig.level = 0.05,insig = "blank",number.cex = 0.3
        ,addrect = 4)
    #dev.off()

  }
    
 

#' create the model using caret package 
#'
#' This function train the model using 80% of the data and retrive the ROC using the 20% of the data
#' @param data which you want to build the model for
#' @param train_index theis index was created using machine_learning function and it is required to make sure 
#' that we train and test the data on the same set
#' @param method which method you want to train the model such as 'lda' and 'gbm'
#' @keywords creat the model
#' @export
#' @examples createthemodel(data,train_index,method)
#' createthemodel(data,train_index,method)
#'
#'
#'    
    
createthemodel <-function(data,train_index,method){
    res=list()
     control <- trainControl(method="cv", number=10,classProbs = TRUE,summaryFunction = twoClassSummary)
     irisTrain <- data[ train_index,]
     irisTest  <- data[-train_index,]
     set.seed(7) 
  garbage <- suppressWarnings(capture.output(fit <- train(subtype~., data=irisTrain, 
                                                          method = method, trControl=control,metric="ROC")) )
  cartClasses <- predict( fit, newdata = irisTest,type="prob")
  cartClasses1 <- predict( fit, newdata = irisTest)
  cartConfusion=confusionMatrix(data = cartClasses1, irisTest$subtype)
  #plot(plot(varImp(fit, scale = FALSE,top=20),main=method))
    
    ROC <- roc(predictor=as.numeric(unlist(cartClasses[1])),response=irisTest$subtype,
               levels=rev(levels(irisTest$subtype)),smooth = FALSE)   
    res$model=fit
    res$ROC=ROC
    return(res)
    
 }