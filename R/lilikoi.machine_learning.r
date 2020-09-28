#' A machine learning Function
#'
#' This function for classification using 8 different machine learning algorithms
#' and it plots the ROC curves and the AUC, SEN, and specificty
#'
#' @param MLmatrix selected pathway deregulation score or metabolites expression matrix
#' @param measurementLabels measurement label for samples
#' @param significantPathways selected pathway names
#' @param trainportion train percentage of the total sample size
#' @param cvnum number of folds
#' @param dlround epoch number for the deep learning method
#' @param Rpart TRUE if run Rpart method
#' @param LDA TRUE if run LDA method
#' @param SVM TRUE if run SVM method
#' @param RF TRUE if run random forest method
#' @param GBM TRUE if run GBM method
#' @param PAM TRUE if run PAM method
#' @param LOG TRUE if run LOG method
#' @param DL TRUE if run deep learning method
#' @keywords classifcation
#' @import caret gbm scales
#' @importFrom graphics legend par
#' @importFrom utils capture.output
#' @importFrom pROC auc roc smooth
#' @importFrom glmnet glmnet
#' @importFrom reshape melt
#' @importFrom Metrics accuracy
#' @importFrom MLmetrics Accuracy AUC Precision Recall Sensitivity Specificity
#' @return Evaluation results and plots of all 8 machine learning algorithms, along with variable importance plots.
#' @export
#' @examples
#' \donttest{
#' lilikoi.machine_learning(MLmatrix = Metadata, measurementLabels = Metadata$Label,
#' significantPathways = 0,
#' trainportion = 0.8, cvnum = 10, dlround=50,Rpart=TRUE,
#' LDA=TRUE,SVM=TRUE,RF=TRUE,GBM=TRUE,PAM=TRUE,LOG=TRUE,DL=TRUE)
#' }

lilikoi.machine_learning <- function (MLmatrix = PDSmatrix, measurementLabels = Label,
                              significantPathways = selected_Pathways_Weka,
                              trainportion = 0.8, cvnum = 10, dlround=50,Rpart=TRUE,
                              LDA=TRUE,SVM=TRUE,RF=TRUE,GBM=TRUE,PAM=TRUE,LOG=TRUE,DL=TRUE) {


  unilabels <- unique(measurementLabels)
  MLmatrix$Label <- as.factor(MLmatrix$Label)
  # if we perform metabolites only machine learning, cancer_df = MLmatrix
  if (significantPathways[1]>0){
    prostate_df <- data.frame(t(MLmatrix[significantPathways,]), Label = measurementLabels, check.names = T)
    colnames(prostate_df)[which(names(prostate_df) == "Label")] <- "subtype"
  } else {
    mPDSmatrix <- t(MLmatrix[-1])
    prostate_df <- data.frame(t(mPDSmatrix), Label = measurementLabels, check.rows = TRUE)
    colnames(prostate_df) <- c(row.names(mPDSmatrix), 'subtype')
  }

  prostate_df$subtype <- as.factor(prostate_df$subtype)

  ##################
  # Quantile normalization, each column corresponds to a sample and each row is a metabolites.
  metadatanorm=normalize.quantiles(t(as.matrix(prostate_df[,-ncol(prostate_df)])))
  #################
  prostate_df[1:(ncol(prostate_df)-1)]=t(metadatanorm)

  colnames(prostate_df) <- make.names(colnames(prostate_df))

  response <- 'subtype'
  predictors <- setdiff(names(prostate_df), response)

  n <- Rpart+LDA+SVM+RF+GBM+PAM+LOG+DL

  if(length(unique(measurementLabels)) == 2){
    performance_training = matrix(rep(0, len = 8*n), nrow = 8)
    performance_testing = matrix(rep(0, len = 8*n), nrow = 8)
    # performance = matrix(rep(0, len = 64), nrow = 8)
  }else{
    performance_training = matrix(rep(0, len = 7*n), nrow = 7)
    performance_testing = matrix(rep(0, len = 7*n), nrow = 7)
    # performance = matrix(rep(0, len = 64), nrow = 7)
  }


  performance_training_list <- list()
  performance_testing_list <- list()
  perfromance_list <- list()
  method_list <- list()
  method <- c()

  # model <- list()

  for (k in 1:10){

    ###############Shuffle stat first
    rand <- sample(nrow(prostate_df))
    prostate_df=prostate_df[rand, ]

    ###############Randomly Split  the data in to training and testing
    trainIndex <- createDataPartition(prostate_df$subtype, p = .8,list = FALSE,times = 1)
    irisTrain <- prostate_df[ trainIndex,]
    irisTest  <- prostate_df[-trainIndex,]
    # irisTrain$subtype=as.factor(paste0("X",irisTrain$subtype))
    # irisTest$subtype=as.factor(paste0("X",irisTest$subtype))

    ###### Create a new summaryFunction called twoClassSummary2, based on twoClassSummary.
    MySummary  <- function(data, lev = NULL, model = NULL){
      a1 <- defaultSummary(data, lev, model)
      b1 <- twoClassSummary(data, lev, model)
      c1 <- prSummary(data, lev, model)
      d1 <- multiClassSummary(data, lev, model)
      out <- c(a1, b1, c1, d1)

      return(out)
    }
    control <- trainControl(method = "cv", number = cvnum, classProbs = TRUE,
                            summaryFunction = MySummary)

    idx <- 0  # Set index for each ML model in order to save results

    #Rpart####
    if (Rpart == TRUE){
      idx <- idx + 1
      method <- c(method, "Rpart")

      set.seed(7)
      # fit.cart <- train(subtype~., data=irisTrain, method = 'rpart', trControl=control,metric="ROC") #loclda
      garbage <- capture.output(fit.cart <- train(subtype ~ .,
                                                  data = irisTrain, method = "rpart", trControl = control,
                                                  metric = "ROC"))
      # model[[idx]] = fit.cart

      performance_training[1,idx]=max(fit.cart$results$ROC)#AUC
      performance_training[2,idx]=fit.cart$results$Sens[which.max(fit.cart$results$ROC)]# sen
      performance_training[3,idx]=fit.cart$results$Spec[which.max(fit.cart$results$ROC)]# spec
      performance_training[4,idx]=fit.cart$results$Accuracy[which.max(fit.cart$results$ROC)]#accuracy
      performance_training[5,idx]=fit.cart$results$Precision[which.max(fit.cart$results$ROC)]#precision
      performance_training[6,idx]=fit.cart$results$Recall[which.max(fit.cart$results$ROC)]#recall = sens
      performance_training[7,idx]=fit.cart$results$F1[which.max(fit.cart$results$ROC)]#F1

      if(length(unilabels) == 2){
        performance_training[8,idx]=fit.cart$results$Balanced_Accuracy[which.max(fit.cart$results$ROC)]#BALANCED ACCURACY
      }


      ##Get the 5 evaluation metrics from testing data
      if(nrow(irisTest) > 0){

        cartClasses <- predict(fit.cart, newdata = irisTest, type = "prob")
        cartClasses1 <- predict(fit.cart, newdata = irisTest)
        cartConfusion = confusionMatrix(data = cartClasses1, irisTest$subtype)
        cart.ROC <- roc(predictor = cartClasses[,1], response = irisTest$subtype,
                        levels = rev(levels(irisTest$subtype)))
        performance_testing[1,idx]=as.numeric(cart.ROC$auc)#AUC
        performance_testing[2,idx]=cartConfusion$byClass[1]#SENS
        performance_testing[3,idx]=cartConfusion$byClass[2]#SPEC
        performance_testing[4,idx]=cartConfusion$overall[1]#accuracy
        performance_testing[5,idx]=cartConfusion$byClass[5]#precision
        performance_testing[6,idx]=cartConfusion$byClass[6]#recall = sens
        performance_testing[7,idx]=cartConfusion$byClass[7]#F1

        if(length(unilabels) == 2){
          performance_testing[8,idx]=cartConfusion$byClass[11]#BALANCED ACCURACY
        }
      }
    }




    #LDA####
    if (LDA == TRUE){
      idx <- idx + 1
      method <- c(method, "LDA")

      set.seed(7)
      fit.lda <- train(subtype~., data=irisTrain, method = 'lda', trControl=control,metric="ROC") #loclda

      # model[[idx]] = fit.lda
      performance_training[1,idx]=max(fit.lda$results$ROC)#AUC
      performance_training[2,idx]=fit.lda$results$Sens[which.max(fit.lda$results$ROC)]# sen
      performance_training[3,idx]=fit.lda$results$Spec[which.max(fit.lda$results$ROC)]# spec
      performance_training[4,idx]=fit.lda$results$Accuracy[which.max(fit.lda$results$ROC)]#accuracy
      performance_training[5,idx]=fit.lda$results$Precision[which.max(fit.lda$results$ROC)]#precision
      performance_training[6,idx]=fit.lda$results$Recall[which.max(fit.lda$results$ROC)]#recall = sens
      performance_training[7,idx]=fit.lda$results$F1[which.max(fit.lda$results$ROC)]#F1
      if(length(unilabels) == 2){
        performance_training[8,idx]=fit.lda$results$Balanced_Accuracy[which.max(fit.lda$results$ROC)]#BALANCED ACCURACY
      }


      ##Get the 5 evaluation metrics from testing data
      if(nrow(irisTest) > 0){

        ldaClasses <- predict( fit.lda, newdata = irisTest,type="prob")
        ldaClasses1 <- predict( fit.lda, newdata = irisTest)
        ldaConfusion=confusionMatrix(data = ldaClasses1, irisTest$subtype)
        lda.ROC <- roc(predictor=ldaClasses[,1],response=irisTest$subtype,levels=rev(levels(irisTest$subtype)))
        performance_testing[1,idx]=as.numeric(lda.ROC$auc)#AUC
        performance_testing[2,idx]=ldaConfusion$byClass[1]#SENS
        performance_testing[3,idx]=ldaConfusion$byClass[2]#SPEC
        performance_testing[4,idx]=ldaConfusion$overall[1]#accuracy
        performance_testing[5,idx]=ldaConfusion$byClass[5]#precision
        performance_testing[6,idx]=ldaConfusion$byClass[6]#recall = sens
        performance_testing[7,idx]=ldaConfusion$byClass[7]#F1
        if(length(unilabels) == 2){
          performance_testing[8,idx]=ldaConfusion$byClass[11]#BALANCED ACCURACY
        }
      }

    }




    #SVM####
    if (SVM == TRUE){
      idx <- idx + 1
      method <- c(method, "SVM")

      set.seed(7)
      garbage <- capture.output(fit.svm <- train(subtype ~ ., data = irisTrain,
                                                 method = "svmRadial", trControl = control, metric = "ROC"))
      # model[[idx]] = fit.svm
      performance_training[1,idx]=max(fit.svm$results$ROC) #AUC
      performance_training[2,idx]=fit.svm$results$Sens[which.max(fit.svm$results$ROC)]# sen
      performance_training[3,idx]=fit.svm$results$Spec[which.max(fit.svm$results$ROC)]# spec
      performance_training[4,idx]=fit.svm$results$Accuracy[which.max(fit.svm$results$ROC)]#accuracy
      performance_training[5,idx]=fit.svm$results$Precision[which.max(fit.svm$results$ROC)]#precision
      performance_training[6,idx]=fit.svm$results$Recall[which.max(fit.svm$results$ROC)]#recall = sens
      performance_training[7,idx]=fit.svm$results$F1[which.max(fit.svm$results$ROC)]#F1
      if(length(unilabels) == 2){
        performance_training[8,idx]=fit.svm$results$Balanced_Accuracy[which.max(fit.svm$results$ROC)]#BALANCED ACCURACY
      }

      if(nrow(irisTest) > 0){

        svmClasses <- predict( fit.svm, newdata = irisTest,type="prob")
        svmClasses1 <- predict( fit.svm, newdata = irisTest)
        svmConfusion=confusionMatrix(data = svmClasses1, irisTest$subtype)
        svm.ROC <- roc(predictor=svmClasses[,1],response=irisTest$subtype,levels=rev(levels(irisTest$subtype)))
        performance_testing[1,idx]=as.numeric(svm.ROC$auc)#AUC
        performance_testing[2,idx]=svmConfusion$byClass[1]#SENS
        performance_testing[3,idx]=svmConfusion$byClass[2]#SPEC
        performance_testing[4,idx]=svmConfusion$overall[1]#accuracy
        performance_testing[5,idx]=svmConfusion$byClass[5]#precision
        performance_testing[6,idx]=svmConfusion$byClass[6]#recall = sens
        performance_testing[7,idx]=svmConfusion$byClass[7]#F1
        if(length(unilabels) == 2){
          performance_testing[8,idx]=svmConfusion$byClass[11]#BALANCED ACCURACY
        }
      }

    }


    #RF####
    if (RF ==TRUE){
      idx <- idx + 1
      method <- c(method, "RF")

      set.seed(7)
      garbage <- capture.output(fit.rf <- train(subtype ~ ., data = irisTrain,
                                                method = "rf", trControl = control, metric = "ROC"))

      # model[[idx]] = fit.rf
      performance_training[1,idx]=max(fit.rf$results$ROC) #AUC
      performance_training[2,idx]=fit.rf$results$Sens[which.max(fit.rf$results$ROC)]# sen
      performance_training[3,idx]=fit.rf$results$Spec[which.max(fit.rf$results$ROC)]# spec
      performance_training[4,idx]=fit.rf$results$Accuracy[which.max(fit.rf$results$ROC)]#accuracy
      performance_training[5,idx]=fit.rf$results$Precision[which.max(fit.rf$results$ROC)]#precision
      performance_training[6,idx]=fit.rf$results$Recall[which.max(fit.rf$results$ROC)]#recall = sens
      performance_training[7,idx]=fit.rf$results$F1[which.max(fit.rf$results$ROC)]#F1
      if(length(unilabels) == 2){
        performance_training[8,idx]=fit.rf$results$Balanced_Accuracy[which.max(fit.rf$results$ROC)]#BALANCED ACCURACY
      }

      # Get testing metrics
      if(nrow(irisTest) > 0){

        rfClasses <- predict( fit.rf, newdata = irisTest,type="prob")
        rfClasses1 <- predict( fit.rf, newdata = irisTest)
        rfConfusion=confusionMatrix(data = rfClasses1, irisTest$subtype)
        rf.ROC <- roc(predictor=rfClasses[,1],response=irisTest$subtype,levels=rev(levels(irisTest$subtype)))
        performance_testing[1,idx]=as.numeric(rf.ROC$auc)#AUC
        performance_testing[2,idx]=rfConfusion$byClass[1]#SENS
        performance_testing[3,idx]=rfConfusion$byClass[2]#SPEC
        performance_testing[4,idx]=rfConfusion$overall[1]#accuracy
        performance_testing[5,idx]=rfConfusion$byClass[5]#precision
        performance_testing[6,idx]=rfConfusion$byClass[6]#recall = sens
        performance_testing[7,idx]=rfConfusion$byClass[7]#F1
        if(length(unilabels) == 2){
          performance_testing[8,idx]=rfConfusion$byClass[11]#BALANCED ACCURACY
        }
      }
    }


    #GBM####
    if (GBM ==TRUE){
      idx <- idx + 1
      method <- c(method, "GBM")

      set.seed(7)
      # fit.gbm <- train(subtype~., data=irisTrain, method="gbm", trControl=control,metric="ROC")
      garbage <- suppressWarnings(capture.output(fit.gbm <- train(subtype ~
                                                                    ., data = irisTrain, method = "gbm", trControl = control,
                                                                  metric = "ROC")))


      # model[[idx]] = fit.gbm
      performance_training[1,idx]=max(fit.gbm$results$ROC) #AUC
      performance_training[2,idx]=fit.gbm$results$Sens[which.max(fit.gbm$results$ROC)]# sen
      performance_training[3,idx]=fit.gbm$results$Spec[which.max(fit.gbm$results$ROC)]# spec
      performance_training[4,idx]=fit.gbm$results$Accuracy[which.max(fit.gbm$results$ROC)]#accuracy
      performance_training[5,idx]=fit.gbm$results$Precision[which.max(fit.gbm$results$ROC)]#precision
      performance_training[6,idx]=fit.gbm$results$Recall[which.max(fit.gbm$results$ROC)]#recall = sens
      performance_training[7,idx]=fit.gbm$results$F1[which.max(fit.gbm$results$ROC)]#F1
      if(length(unilabels) == 2){
        performance_training[8,idx]=fit.gbm$results$Balanced_Accuracy[which.max(fit.gbm$results$ROC)]#BALANCED ACCURACY
      }

      if(nrow(irisTest) > 0){

        gbmClasses <- predict( fit.gbm, newdata = irisTest,type="prob")
        gbmClasses1 <- predict( fit.gbm, newdata = irisTest)
        gbmConfusion=confusionMatrix(data = gbmClasses1, irisTest$subtype)
        gbm.ROC <- roc(predictor=gbmClasses[,1],response=irisTest$subtype,levels=rev(levels(irisTest$subtype)))
        performance_testing[1,idx]=as.numeric(gbm.ROC$auc)#AUC
        performance_testing[2,idx]=gbmConfusion$byClass[1]#SENS
        performance_testing[3,idx]=gbmConfusion$byClass[2]#SPEC
        performance_testing[4,idx]=gbmConfusion$overall[1]#accuracy
        performance_testing[5,idx]=gbmConfusion$byClass[5]#precision
        performance_testing[6,idx]=gbmConfusion$byClass[6]#recall = sens
        performance_testing[7,idx]=gbmConfusion$byClass[7]#F1
        if(length(unilabels) == 2){
          performance_testing[8,idx]=gbmConfusion$byClass[11]#BALANCED ACCURACY
        }
      }
    }


    #PAM####
    if (PAM == TRUE){
      idx <- idx + 1
      method <- c(method, "PAM")

      set.seed(7)
      fit.pam <- train(subtype~., data=irisTrain, method="pam", trControl=control,metric="ROC")#plr

      # model[[idx]] = fit.pam
      performance_training[1,idx]=max(fit.pam$results$ROC) #AUC
      performance_training[2,idx]=fit.pam$results$Sens[which.max(fit.pam$results$ROC)]# sen
      performance_training[3,idx]=fit.pam$results$Spec[which.max(fit.pam$results$ROC)]# spec
      performance_training[4,idx]=fit.pam$results$Accuracy[which.max(fit.pam$results$ROC)]#accuracy
      performance_training[5,idx]=fit.pam$results$Precision[which.max(fit.pam$results$ROC)]#precision
      performance_training[6,idx]=fit.pam$results$Recall[which.max(fit.pam$results$ROC)]#recall = sens
      performance_training[7,idx]=fit.pam$results$F1[which.max(fit.pam$results$ROC)]#F1
      if(length(unilabels) == 2){
        performance_training[8,idx]=fit.pam$results$Balanced_Accuracy[which.max(fit.pam$results$ROC)]#BALANCED ACCURACY
      }


      if(nrow(irisTest) > 0){

        pamClasses <- predict(fit.pam, newdata = irisTest,type="prob")
        pamClasses1 <- predict(fit.pam, newdata = irisTest)
        pamConfusion=confusionMatrix(data = pamClasses1, irisTest$subtype)
        pam.ROC <- roc(predictor=pamClasses[,1],response=irisTest$subtype,levels=rev(levels(irisTest$subtype)))
        performance_testing[1,idx]=as.numeric(pam.ROC$auc)#AUC
        performance_testing[2,idx]=pamConfusion$byClass[1]#SENS
        performance_testing[3,idx]=pamConfusion$byClass[2]#SPEC
        performance_testing[4,idx]=pamConfusion$overall[1]#accuracy
        performance_testing[5,idx]=pamConfusion$byClass[5]#precision
        performance_testing[6,idx]=pamConfusion$byClass[6]#recall = sens
        performance_testing[7,idx]=pamConfusion$byClass[7]#F1
        if(length(unilabels) == 2){
          performance_testing[8,idx]=pamConfusion$byClass[11]#BALANCED ACCURACY
        }
      }
    }


    #LOG####
    if (LOG == TRUE){
      idx <- idx + 1
      method <- c(method, "LOG")

      set.seed(7)
      fit.log <- train(subtype~., data=irisTrain, method="glmnet", trControl=control,metric="ROC")#plr

      # model[[idx]] = fit.log
      performance_training[1,idx]=max(fit.log$results$ROC) #AUC
      performance_training[2,idx]=fit.log$results$Sens[which.max(fit.log$results$ROC)]# sen
      performance_training[3,idx]=fit.log$results$Spec[which.max(fit.log$results$ROC)]# spec
      performance_training[4,idx]=fit.log$results$Accuracy[which.max(fit.log$results$ROC)]#accuracy
      performance_training[5,idx]=fit.log$results$Precision[which.max(fit.log$results$ROC)]#precision
      performance_training[6,idx]=fit.log$results$Recall[which.max(fit.log$results$ROC)]#recall = sens
      performance_training[7,idx]=fit.log$results$F1[which.max(fit.log$results$ROC)]#F1
      if(length(unilabels) == 2){
        performance_training[8,idx]=fit.log$results$Balanced_Accuracy[which.max(fit.log$results$ROC)]#BALANCED ACCURACY
      }

      # Get test metrics
      if(nrow(irisTest) > 0){

        logClasses <- predict( fit.log, newdata = irisTest,type="prob")
        logClasses1 <- predict( fit.log, newdata = irisTest)
        logConfusion=confusionMatrix(data = logClasses1, irisTest$subtype)
        log.ROC <- roc(predictor=logClasses[,1],response=irisTest$subtype,levels=rev(levels(irisTest$subtype)))
        performance_testing[1,idx]=as.numeric(log.ROC$auc)#AUC
        performance_testing[2,idx]=logConfusion$byClass[1]#SENS
        performance_testing[3,idx]=logConfusion$byClass[2]#SPEC
        performance_testing[4,idx]=logConfusion$overall[1]#accuracy
        performance_testing[5,idx]=logConfusion$byClass[5]#precision
        performance_testing[6,idx]=logConfusion$byClass[6]#recall = sens
        performance_testing[7,idx]=logConfusion$byClass[7]#F1
        if(length(unilabels) == 2){
          performance_testing[8,idx]=logConfusion$byClass[11]#BALANCED ACCURACY
        }
      }
    }



    #Deep learning#######
    if (DL == TRUE){
      idx <- idx + 1
      method <- c(method, "DL")

      localH2O = h2o.init()
      irisTrain$subtype <- as.factor(irisTrain$subtype)
      prostate.hex<-as.h2o(irisTrain, destination_frame="train.hex")
      #valid <- as.h2o(irisTest, destination_frame="test.hex")
      #prostate.hex$subtype <- as.factor(prostate.hex$subtype) ##make categorical
      hyper_params <- list(
        activation=c("Rectifier","Tanh"),
        hidden=list(c(100),c(200),c(10,10),c(20,20),c(50,50),c(30,30,30),c(25,25,25,25)),
        input_dropout_ratio=c(0,0.05,0.1),
        #hidden_dropout_ratios=c(0.6,0.5,0.6,0.6),
        l1=seq(0,1e-4,1e-6),
        l2=seq(0,1e-4,1e-6),
        train_samples_per_iteration =c(0,-2),
        epochs = c(dlround),
        momentum_start=c(0,0.5),
        rho=c(0.5,0.99),
        quantile_alpha=c(0,1),
        huber_alpha=seq(0,1) )

      search_criteria = list(strategy = "RandomDiscrete", max_models = 100, stopping_rounds=5,
                             stopping_tolerance=1e-2)
      dl_random_grid <- h2o.grid(
        algorithm="deeplearning",
        #grid_id = paste("dl_do3",k,sep = "_"),
        grid_id = "dl_grid_randome1",
        training_frame=prostate.hex,
        #validation_frame=valid,
        #hidden=c(25,25,25,25),
        x=predictors,
        y="subtype",
        #pretrained_autoencoder="dl_grid_random1_model_42",
        #autoencoder = TRUE,
        seed=7,
        #adaptive_rate=T, ## manually tuned learning rate
        #momentum_start=0.5, ## manually tuned momentum
        #momentum_stable=0.9,
        #momentum_ramp=1e7,
        variable_importances=TRUE,
        export_weights_and_biases=T,
        standardize=T,
        stopping_metric="misclassification",
        stopping_tolerance=1e-2, ## stop when logloss does not improve by >=1% for 2 scoring events
        stopping_rounds=2,
        #score_validation_samples=10000, ## downsample validation set for faster scoring
        score_duty_cycle=0.025, ## don't score more than 2.5% of the wall time
        #max_w2=10, ## can help improve stability for Rectifier
        hyper_params = hyper_params,
        search_criteria = search_criteria,
        nfolds=10
      )

      grid <- h2o.getGrid("dl_grid_randome1",sort_by="mse",decreasing=FALSE)
      grid@summary_table[1,]
      best_model <- h2o.getModel(grid@model_ids[[1]]) ## model with lowest logloss
      # model[[idx]] <- best_model
      performance_training[1,idx]=as.numeric(best_model@model$cross_validation_metrics_summary$mean)[2] #AUC
      performance_training[2,idx]=as.numeric(best_model@model$cross_validation_metrics_summary$mean)[17]# sen
      performance_training[3,idx]=as.numeric(best_model@model$cross_validation_metrics_summary$mean)[19]# spec

      performance_training[4,idx]=as.numeric(best_model@model$cross_validation_metrics_summary$mean)[1] #accuracy
      performance_training[5,idx]=as.numeric(best_model@model$cross_validation_metrics_summary$mean)[15]#precision
      performance_training[6,idx]=as.numeric(best_model@model$cross_validation_metrics_summary$mean)[17]#recall
      performance_training[7,idx]=as.numeric(best_model@model$cross_validation_metrics_summary$mean)[6]# f1
      if(length(unilabels) == 2){
        performance_training[8,idx]=(performance_training[2,idx]+performance_training[3,idx])/2
      }

      # Get test metrics
      if(nrow(irisTest) > 0){
        # irisTest$subtype <- as.factor(irisTest$subtype)
        perf=h2o.performance(best_model,as.h2o(irisTest, destination_frame="test.hex"))
        performance_testing[1,idx]=as.numeric(h2o.auc(perf,h2o.find_threshold_by_max_metric(perf,"f1"))[[1]])#AUC
        performance_testing[2,idx]=as.numeric(h2o.sensitivity(perf,h2o.find_threshold_by_max_metric(perf,"f1"))[[1]])#SENS
        performance_testing[3,idx]=as.numeric(h2o.specificity(perf,h2o.find_threshold_by_max_metric(perf,"f1"))[[1]])#SPEC
        performance_testing[4,idx]=as.numeric(h2o.accuracy(perf,h2o.find_threshold_by_max_metric(perf,"f1"))[[1]])#accuracy
        performance_testing[5,idx]=as.numeric(h2o.precision(perf,h2o.find_threshold_by_max_metric(perf,"f1")[[1]]))#precision
        performance_testing[6,idx]=as.numeric(h2o.sensitivity(perf,h2o.find_threshold_by_max_metric(perf,"f1"))[[1]])#recall = sens
        performance_testing[7,idx]=as.numeric(h2o.F1(perf,h2o.find_threshold_by_max_metric(perf,"f1"))[[1]])#F1

        if(length(unilabels) == 2){
          performance_testing[8,idx]=(performance_testing[2,idx]+performance_testing[3,idx])/2 #BALANCED ACCURACY
        }
      }

    }

      performance_testing_list[[k]] <- performance_testing
      performance_training_list[[k]] <- performance_training
      method_list[[k]] <- method


      if(length(unique(measurementLabels)) == 2){
        performance_training = matrix(rep(0, len = 8*n), nrow = 8)
        performance_testing = matrix(rep(0, len = 8*n), nrow = 8)
      }else{
        performance_training = matrix(rep(0, len = 7*n), nrow = 7)
        performance_testing = matrix(rep(0, len = 7*n), nrow = 7)
      }
      method = c()
    }

    list_test <- performance_testing_list
    list_train <- performance_training_list
    list_method <- method_list[[1]]


    AUC_train <- lapply(list_train, function(x) x[1, ])
    AUC_test <- lapply(list_test, function(x) x[1, ])

    SENS_train <- lapply(list_train, function(x) x[2, ])
    SENS_test <- lapply(list_test, function(x) x[2, ])

    SPEC_train <- lapply(list_train, function(x) x[3, ])
    SPEC_test <- lapply(list_test, function(x) x[3, ])

    # accuracy_test <- lapply(list_test, function(x) x[4, ])
    # precision_test <- lapply(list_test, function(x) x[5, ])
    # recall_test <- lapply(list_test, function(x) x[6, ])
    F1_train <- lapply(list_train, function(x) x[7, ])
    F1_test <- lapply(list_test, function(x) x[7, ])

    Balanced_accuracy_train <- lapply(list_train, function(x) x[8, ])
    Balanced_accuracy_test <- lapply(list_test, function(x) x[8, ])

    output1 <- do.call(rbind, lapply(AUC_train, matrix, ncol = n,
                                     byrow = TRUE))
    output2 <- do.call(rbind, lapply(AUC_test, matrix, ncol = n,
                                     byrow = TRUE))
    output3 <- do.call(rbind, lapply(SENS_train, matrix, ncol = n,
                                     byrow = TRUE))
    output4 <- do.call(rbind, lapply(SENS_test, matrix, ncol = n,
                                     byrow = TRUE))
    output5 <- do.call(rbind, lapply(SPEC_train, matrix, ncol = n,
                                     byrow = TRUE))
    output6 <- do.call(rbind, lapply(SPEC_test, matrix, ncol = n,
                                     byrow = TRUE))
    output7 <- do.call(rbind, lapply(F1_train, matrix, ncol = n,
                                     byrow = TRUE))
    output8 <- do.call(rbind, lapply(F1_test, matrix, ncol = n,
                                     byrow = TRUE))
    output9 <- do.call(rbind, lapply(Balanced_accuracy_train,
                                     matrix, ncol = n, byrow = TRUE))
    output10 <- do.call(rbind, lapply(Balanced_accuracy_test,
                                     matrix, ncol = n, byrow = TRUE))

    AUC_train_mean <- apply(output1, 2, mean)
    AUC_train_sd <- apply(output1, 2, sd)
    AUC_test_mean <- apply(output2, 2, mean)
    AUC_test_sd <- apply(output2, 2, sd)

    SENS_train_mean <- apply(output3, 2, mean)
    SENS_train_sd <- apply(output3, 2, sd)
    SENS_test_mean <- apply(output4, 2, mean)
    SENS_test_sd <- apply(output4, 2, sd)

    SPEC_train_mean=apply(output5,2,mean)
    SPEC_train_sd=apply(output5,2,sd)
    SPEC_test_mean=apply(output6,2,mean)
    SPEC_test_sd=apply(output6,2,sd)

    F1_train_mean <- apply(output7, 2, mean)
    F1_train_sd <- apply(output7, 2, sd)
    F1_test_mean <- apply(output8, 2, mean)
    F1_test_sd <- apply(output8, 2, sd)

    Balanced_accuracy_train_mean <- apply(output9, 2, mean)
    Balanced_accuracy_train_sd <- apply(output9, 2, sd)
    Balanced_accuracy_test_mean <- apply(output10, 2, mean)
    Balanced_accuracy_test_sd <- apply(output10, 2, sd)

    trainingORtesting <- t(cbind(t(rep("training", n)), t(rep("testing", n))))
    testing_only <- t(t(rep("testing", n)))
    training_only <- t(t(rep('training', n)))

    performance_data_test <- data.frame(Algorithm = (rep(t(list_method), 1)),
                                        testing_only,
                                        AUC = data.frame(AUC = t((t(AUC_test_mean)))),
                                        SENS = data.frame(SENS = t((t(SENS_test_mean)))),
                                        SPEC = data.frame(SPEC = t((t(SPEC_test_mean)))),
                                        F1 = data.frame(F1 = t((t(F1_test_mean)))),
                                        Balanced_accuracy = data.frame(Balanced_accuracy = t((t(Balanced_accuracy_test_mean)))))

    performance_data_train <- data.frame(Algorithm = (rep(t(list_method), 1)),
                                         training_only,
                                         AUC = data.frame(AUC = t((t(AUC_train_mean)))),
                                         SENS = data.frame(SENS = t((t(SENS_train_mean)))),
                                         SPEC = data.frame(SPEC = t((t(SPEC_train_mean)))),
                                         F1 = data.frame(F1 = t((t(F1_train_mean)))),
                                         Balanced_accuracy = data.frame(Balanced_accuracy = t((t(Balanced_accuracy_train_mean)))))


    # library(reshape)
    melted_performance_data_test <- suppressMessages(melt(performance_data_test))
    melted_performance_data_train <- suppressMessages(melt(performance_data_train))

    sd_train <- data.frame(sd = rbind(t(t(AUC_train_sd)),t(t(SENS_train_sd)),t(t(SPEC_train_sd)),
                                      t(t(F1_train_sd)),t(t(Balanced_accuracy_train_sd))))
    sd_test <- data.frame(sd = rbind(t(t(AUC_test_sd)),t(t(SENS_test_sd)),t(t(SPEC_test_sd)),
                                     t(t(F1_test_sd)),t(t(Balanced_accuracy_test_sd))))
    melted_performance_data_train$sd <- sd_train
    melted_performance_data_test$sd <- sd_test

    melted_performance_data_train$sd <- unlist(melted_performance_data_train$sd)
    melted_performance_data_test$sd <- unlist(melted_performance_data_test$sd)

    p1 <- ggplot(data = melted_performance_data_train,
                 aes(x = Algorithm, y = value, fill = variable)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                    position=position_dodge(.9)) +
      xlab("") + ylab("") + ggtitle("Training") +
      theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 15, face = "bold"),
            axis.title = element_text(size = 14, face = "bold")) + labs(fill = "")
    # print(p1)

    if(nrow(irisTest) > 0){

      p2 <- ggplot(data = melted_performance_data_test,
                   aes(x = Algorithm, y = value, fill = variable)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                      position=position_dodge(.9)) +
        xlab("") + ylab("") + ggtitle("Testing") +
        theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 15, face = "bold"),
              axis.title = element_text(size = 14, face = "bold")) + labs(fill = "")
      # print(p2)

    }




    plots <- list()
    performance <- list()

    performance$TrainingPerformance <- performance_data_train
    performance$TestingPerformance <- performance_data_test

    plots$TrainingPlot <- p1
    plots$TestingPlot <- p2


    return(c(performance, plots))


  }

