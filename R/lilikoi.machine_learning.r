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
#' @keywords classifcation
#' @import caret gbm scales
#' @importFrom graphics legend par
#' @importFrom utils capture.output
#' @importFrom h2o h2o.init as.h2o h2o.grid h2o.getGrid h2o.getModel h2o.performance h2o.accuracy h2o.precision h2o.sensitivity h2o.F1 h2o.find_threshold_by_max_metric h2o.auc h2o.varimp_plot h2o.deeplearning h2o.getFrame h2o.predict
#' @importFrom pROC auc roc smooth
#' @importFrom glmnet glmnet
#' @importFrom reshape melt
#' @importFrom Metrics accuracy
#' @return Evaluation results and plots of all 8 machine learning algorithms, along with variable importance plots.
#' @export
#' @examples
#' \donttest{
#'  measurementLabels = Metadata$Label
#'  MLmatrix = Metadata
#'  significantPathways = 0
#'  lilikoi.machine_learning(tPDSmatrix, measurementLabels,significantPathways,trainportion = 0.8,
#'     cvnum = 10)
#' }
#'
#'

lilikoi.machine_learning <- function (MLmatrix = PDSmatrix, measurementLabels = Label,
                              significantPathways = selected_Pathways_Weka,
                              trainportion = 0.8, cvnum = 10) {

  unilabels <- unique(measurementLabels)

  # if we perform metabolites only machine learning, cancer_df = MLmatrix
  if (significantPathways[1]>0){
    prostate_df <- data.frame(t(MLmatrix[significantPathways,]), Label = measurementLabels, check.names = T)
  } else {
    prostate_df <- MLmatrix
  }

  colnames(prostate_df)[which(names(prostate_df) == "Label")] <- "subtype"

  ##################
  # Quantile normalization, each column corresponds to a sample and each row is a metabolites.
  metadatanorm=normalize.quantiles(t(as.matrix(prostate_df[,-1])))
  #################
  prostate_df[2:ncol(prostate_df)]=t(metadatanorm)
  response <- 'subtype'
  predictors <- setdiff(names(prostate_df), response)

  prostate_df$subtype <- as.factor(prostate_df$subtype)
  colnames(prostate_df) <- make.names(colnames(prostate_df))

  if(length(unique(measurementLabels)) == 2){
    performance_training = matrix(rep(0, len = 24), nrow = 3)
    performance_testing = matrix(rep(0, len = 64), nrow = 8)
    # performance = matrix(rep(0, len = 64), nrow = 8)
  }else{
    performance_training = matrix(rep(0, len = 21), nrow = 3)
    performance_testing = matrix(rep(0, len = 64), nrow = 7)
    # performance = matrix(rep(0, len = 64), nrow = 7)
  }


  performance_training_list <- list()
  performance_testing_list <- list()
  perfromance_list <- list()

  model <- list()

  for (k in 1:1){
    set.seed(2000)

    ###############Shuffle stat first
    rand <- sample(nrow(prostate_df))
    prostate_df=prostate_df[rand, ]

    ###############Randomly Split  the data in to training and testing
    set.seed(2000)
    trainIndex <- createDataPartition(prostate_df$subtype, p = .8,list = FALSE,times = 1)
    irisTrain <- prostate_df[ trainIndex,]
    irisTest  <- prostate_df[-trainIndex,]


    ###### Create a new summaryFunction called twoClassSummary2, based on twoClassSummary.
    # twoClassSummary2 <- function (data, lev = NULL, model = NULL)
    # {
    #   # if (length(lev) > 2) {
    #   #   stop(paste("Your outcome has", length(lev), "levels. The twoClassSummary() function isn't appropriate."))
    #   # }
    #   # requireNamespaceQuietStop("pROC")
    #   # if (!all(levels(data[, "pred"]) == lev)) {
    #   #   stop("levels of observed and predicted data do not match")
    #   # }
    #   rocObject <- try(pROC::roc(data$obs, data[, lev[1]], direction = ">",
    #                              quiet = TRUE), silent = TRUE)
    #   rocAUC <- if (inherits(rocObject, "try-error"))
    #     NA
    #   else rocObject$auc
    #
    #   sensitivity <- sensitivity(data[, "pred"], data[, "obs"], lev[1])
    #   specificity <- specificity(data[, "pred"], data[, "obs"], lev[2])
    #   precision <- posPredValue(data$pred, data$obs, positive = "ERp")
    #   recall  <- sensitivity(data$pred, data$obs, postive = "ERp")
    #   accuracy <- accuracy(data$pred, data$obs)
    #   f1_val <- (2 * precision * recall) / (precision + recall)
    #   balanced_accu <- (sensitivity+specificity)/2
    #
    #   out <- c(rocAUC, sensitivity, specificity, accuracy,
    #            precision, recall, f1_val,balanced_accu)
    #   names(out) <- c("ROC", "Sens", "Spec", "Accu","Prec", "Recall","F1","Balanced_Accu")
    #   out
    # }

    control <- trainControl(method = "cv", number = cvnum, classProbs = TRUE,
                            summaryFunction = twoClassSummary)
    # savePredictions = returncvtype)

    #C4.5 Rpart###

    set.seed(7)
    garbage <- capture.output(fit.cart <- train(subtype ~ .,
                                                data = irisTrain, method = "rpart", trControl = control,
                                                metric = "ROC"))
    model[[1]] = fit.cart

    performance_training[1, 1] = max(fit.cart$results$ROC) # AUC
    performance_training[2, 1] = fit.cart$results$Sens[which.max(fit.cart$results$ROC)] # Sensitivity, Recall
    performance_training[3, 1] = fit.cart$results$Spec[which.max(fit.cart$results$ROC)] # Specificity

    # performance_training[4, 1] <- fit.cart$results$Accu[which.max(fit.cart$results$ROC)]
    # performance_training[5, 1] <- fit.cart$results$Prec[which.max(fit.cart$results$ROC)]
    # performance_training[6, 1] <- fit.cart$results$Recall[which.max(fit.cart$results$ROC)]
    # performance_training[7, 1] <- fit.cart$results$F1[which.max(fit.cart$results$ROC)]
    # if(length(unilabels) == 2){
    #   performance_training[8, 1] <- fit.cart$results$Balanced_Accu[which.max(fit.cart$results$ROC)]
    # }

    ##Get the 5 evaluation metrics from testing data
    if(nrow(irisTest) > 0){

      cartClasses <- predict(fit.cart, newdata = irisTest, type = "prob")
      cartClasses1 <- predict(fit.cart, newdata = irisTest)
      cartConfusion = confusionMatrix(data = cartClasses1, irisTest$subtype)
      cart.ROC <- roc(predictor = cartClasses[,1], response = irisTest$subtype,
                      levels = rev(levels(irisTest$subtype)))

      performance_testing[1, 1] = as.numeric(cart.ROC$auc)
      performance_testing[2, 1] = cartConfusion$byClass[1]
      performance_testing[3, 1] = cartConfusion$byClass[2]
      performance_testing[4, 1] = cartConfusion$overall[1]
      performance_testing[5, 1] = cartConfusion$byClass[5]
      performance_testing[6, 1] = cartConfusion$byClass[6]
      performance_testing[7, 1] = cartConfusion$byClass[7]

      if(length(unilabels) == 2){
        performance_testing[8, 1] <- cartConfusion$byClass[11]
      }


    }



    #LDA###
    set.seed(7)
    garbage <- suppressWarnings(capture.output(fit.lda <- train(subtype ~
                                                                  ., data = irisTrain, method = "lda", trControl = control,
                                                                metric = "ROC", trace = F)))
    model[[2]] = fit.lda
    performance_training[1, 2] = max(fit.lda$results$ROC)
    performance_training[2, 2] = fit.lda$results$Sens[which.max(fit.lda$results$ROC)]
    performance_training[3, 2] = fit.lda$results$Spec[which.max(fit.lda$results$ROC)]


    # performance_training[4, 2] <- fit.lda$results$Accu[which.max(fit.lda$results$ROC)]
    # performance_training[5, 2] <- fit.lda$results$Prec[which.max(fit.lda$results$ROC)]
    # performance_training[6, 2] <- fit.lda$results$Recall[which.max(fit.lda$results$ROC)]
    # performance_training[7, 2] <- fit.lda$results$F1[which.max(fit.lda$results$ROC)]
    # if(length(unilabels) == 2){
    #   performance_training[8, 2] <- fit.lda$results$Balanced_Accu[which.max(fit.lda$results$ROC)]
    # }


    ##Get the 5 evaluation metrics from testing data
    if(nrow(irisTest) > 0){

      ldaClasses <- predict(fit.lda, newdata = irisTest, type = "prob")
      ldaClasses1 <- predict(fit.lda, newdata = irisTest)
      ldaConfusion = confusionMatrix(data = ldaClasses1, irisTest$subtype)
      lda.ROC <- roc(predictor = ldaClasses[,1], response = irisTest$subtype,
                     levels = rev(levels(irisTest$subtype)))
      performance_testing[1, 2] = as.numeric(lda.ROC$auc)
      performance_testing[2, 2] = ldaConfusion$byClass[1]
      performance_testing[3, 2] = ldaConfusion$byClass[2]
      performance_testing[4, 2] = ldaConfusion$overall[1]
      performance_testing[5, 2] = ldaConfusion$byClass[5]
      performance_testing[6, 2] = ldaConfusion$byClass[6]
      performance_testing[7, 2] = ldaConfusion$byClass[7]
      if(length(unilabels) == 2){
        performance_testing[8, 2] <- ldaConfusion$byClass[11]
      }


    }

    #SVM###
    set.seed(7)
    garbage <- capture.output(fit.svm <- train(subtype ~ ., data = irisTrain,
                                               method = "svmRadial", trControl = control, metric = "ROC"))
    model[[3]] = fit.svm
    performance_training[1, 3] = max(fit.svm$results$ROC)
    performance_training[2, 3] = fit.svm$results$Sens[which.max(fit.svm$results$ROC)]
    performance_training[3, 3] = fit.svm$results$Spec[which.max(fit.svm$results$ROC)]

    # performance_training[4, 3] <- fit.svm$results$Accu[which.max(fit.svm$results$ROC)]
    # performance_training[5, 3] <- fit.svm$results$Prec[which.max(fit.svm$results$ROC)]
    # performance_training[6, 3] <- fit.svm$results$Recall[which.max(fit.svm$results$ROC)]
    # performance_training[7, 3] <- fit.svm$results$F1[which.max(fit.svm$results$ROC)]
    # if(length(unilabels) == 2){
    #   performance_training[8, 3] <- fit.svm$results$Balanced_Accu[which.max(fit.svm$results$ROC)]
    # }



    if(nrow(irisTest) > 0){

      svmClasses <- predict(fit.svm, newdata = irisTest, type = "prob")
      svmClasses1 <- predict(fit.svm, newdata = irisTest)
      svmConfusion = confusionMatrix(data = svmClasses1, irisTest$subtype)
      svm.ROC <- roc(predictor = svmClasses[,1], response = irisTest$subtype,
                     levels = rev(levels(irisTest$subtype)))
      performance_testing[1, 3] = as.numeric(svm.ROC$auc)
      performance_testing[2, 3] = svmConfusion$byClass[1]
      performance_testing[3, 3] = svmConfusion$byClass[2]
      performance_testing[4, 3] = svmConfusion$overall[1]
      performance_testing[5, 3] = svmConfusion$byClass[5]
      performance_testing[6, 3] = svmConfusion$byClass[6]
      performance_testing[7, 3] = svmConfusion$byClass[7]
      if(length(unilabels) == 2){
        performance_testing[8, 3] <- svmConfusion$byClass[11]
      }


    }

    #RF###
    set.seed(7)
    garbage <- capture.output(fit.rf <- train(subtype ~ ., data = irisTrain,
                                              method = "rf", trControl = control, metric = "ROC"))
    model[[4]] = fit.rf
    performance_training[1, 4] = max(fit.rf$results$ROC)
    performance_training[2, 4] = fit.rf$results$Sens[which.max(fit.rf$results$ROC)]
    performance_training[3, 4] = fit.rf$results$Spec[which.max(fit.rf$results$ROC)]

    # performance_training[4, 4] <- fit.rf$results$Accu[which.max(fit.rf$results$ROC)]
    # performance_training[5, 4] <- fit.rf$results$Prec[which.max(fit.rf$results$ROC)]
    # performance_training[6, 4] <- fit.rf$results$Recall[which.max(fit.rf$results$ROC)]
    # performance_training[7, 4] <- fit.rf$results$F1[which.max(fit.rf$results$ROC)]
    # if(length(unilabels) == 2){
    #   performance_training[8, 4] <- fit.rf$results$Balanced_Accu[which.max(fit.rf$results$ROC)]
    # }



    if(nrow(irisTest) > 0){

      rfClasses <- predict(fit.rf, newdata = irisTest, type = "prob")
      rfClasses1 <- predict(fit.rf, newdata = irisTest)
      rfConfusion = confusionMatrix(data = rfClasses1, irisTest$subtype)
      rf.ROC <- roc(predictor = rfClasses[,1], response = irisTest$subtype,
                    levels = rev(levels(irisTest$subtype)))
      performance_testing[1, 4] = as.numeric(rf.ROC$auc)
      performance_testing[2, 4] = rfConfusion$byClass[1]
      performance_testing[3, 4] = rfConfusion$byClass[2]
      performance_testing[4, 4] = rfConfusion$overall[1]
      performance_testing[5, 4] = rfConfusion$byClass[5]
      performance_testing[6, 4] = rfConfusion$byClass[6]
      performance_testing[7, 4] = rfConfusion$byClass[7]
      if(length(unilabels) == 2){
        performance_testing[8, 4] <- rfConfusion$byClass[11]
      }


    }

    #GBM###
    set.seed(7)
    garbage <- suppressWarnings(capture.output(fit.gbm <- train(subtype ~
                                                                  ., data = irisTrain, method = "gbm", trControl = control,
                                                                metric = "ROC")))

    model[[5]] = fit.gbm
    performance_training[1, 5] = max(fit.gbm$results$ROC)
    performance_training[2, 5] = fit.gbm$results$Sens[which.max(fit.gbm$results$ROC)]
    performance_training[3, 5] = fit.gbm$results$Spec[which.max(fit.gbm$results$ROC)]

    # performance_training[4, 5] <- fit.gbm$results$Accu[which.max(fit.gbm$results$ROC)]
    # performance_training[5, 5] <- fit.gbm$results$Prec[which.max(fit.gbm$results$ROC)]
    # performance_training[6, 5] <- fit.gbm$results$Recall[which.max(fit.gbm$results$ROC)]
    # performance_training[7, 5] <- fit.gbm$results$F1[which.max(fit.gbm$results$ROC)]
    # if(length(unilabels) == 2){
    #   performance_training[8, 5] <- fit.gbm$results$Balanced_Accu[which.max(fit.gbm$results$ROC)]
    # }



    if(nrow(irisTest) > 0){

      gbmClasses <- predict(fit.gbm, newdata = irisTest, type = "prob")
      gbmClasses1 <- predict(fit.gbm, newdata = irisTest)
      gbmConfusion = confusionMatrix(data = gbmClasses1, irisTest$subtype)
      gbm.ROC <- roc(predictor = gbmClasses[,1], response = irisTest$subtype,
                     levels = rev(levels(irisTest$subtype)))
      performance_testing[1, 5] = as.numeric(gbm.ROC$auc)
      performance_testing[2, 5] = gbmConfusion$byClass[1]
      performance_testing[3, 5] = gbmConfusion$byClass[2]
      performance_testing[4, 5] = gbmConfusion$overall[1]
      performance_testing[5, 5] = gbmConfusion$byClass[5]
      performance_testing[6, 5] = gbmConfusion$byClass[6]
      performance_testing[7, 5] = gbmConfusion$byClass[7]
      if(length(unilabels) == 2){
        performance_testing[8, 5] <- gbmConfusion$byClass[11]
      }


    }

    #PAM###
    set.seed(7)
    garbage <- capture.output(fit.pam <- train(subtype ~ ., data = irisTrain,
                                               method = "pam", trControl = control, metric = "ROC"))
    model[[6]] = fit.pam
    performance_training[1, 6] = max(fit.pam$results$ROC)
    performance_training[2, 6] = fit.pam$results$Sens[which.max(fit.pam$results$ROC)]
    performance_training[3, 6] = fit.pam$results$Spec[which.max(fit.pam$results$ROC)]

    # performance_training[4, 6] <- fit.pam$results$Accu[which.max(fit.pam$results$ROC)]
    # performance_training[5, 6] <- fit.pam$results$Prec[which.max(fit.pam$results$ROC)]
    # performance_training[6, 6] <- fit.pam$results$Recall[which.max(fit.pam$results$ROC)]
    # performance_training[7, 6] <- fit.pam$results$F1[which.max(fit.pam$results$ROC)]
    # if(length(unilabels) == 2){
    #   performance_training[8, 6] <- fit.pam$results$Balanced_Accu[which.max(fit.pam$results$ROC)]
    # }



    if(nrow(irisTest) > 0){

      pamClasses <- predict(fit.pam, newdata = irisTest, type = "prob")
      pamClasses1 <- predict(fit.pam, newdata = irisTest)
      pamConfusion = confusionMatrix(data = pamClasses1, irisTest$subtype)
      pam.ROC <- roc(predictor = pamClasses[,1], response = irisTest$subtype,
                     levels = rev(levels(irisTest$subtype)))
      performance_testing[1, 6] = as.numeric(pam.ROC$auc)
      performance_testing[2, 6] = pamConfusion$byClass[1]
      performance_testing[3, 6] = pamConfusion$byClass[2]
      performance_testing[4, 6] = pamConfusion$overall[1]
      performance_testing[5, 6] = pamConfusion$byClass[5]
      performance_testing[6, 6] = pamConfusion$byClass[6]
      performance_testing[7, 6] = pamConfusion$byClass[7]
      if(length(unilabels) == 2){
        performance_testing[8, 6] <- pamConfusion$byClass[11]
      }


    }

    #LOG###
    set.seed(7)
    garbage <- suppressWarnings(capture.output(fit.log <- train(subtype ~
                                                                  ., data = irisTrain, method = "glmnet", trControl = control,
                                                                metric = "ROC")))
    model[[7]] = fit.log
    performance_training[1, 7] = max(fit.log$results$ROC)
    performance_training[2, 7] = fit.log$results$Sens[which.max(fit.log$results$ROC)]
    performance_training[3, 7] = fit.log$results$Spec[which.max(fit.log$results$ROC)]


    # performance_training[4, 7] <- fit.log$results$Accu[which.max(fit.log$results$ROC)]
    # performance_training[5, 7] <- fit.log$results$Prec[which.max(fit.log$results$ROC)]
    # performance_training[6, 7] <- fit.log$results$Recall[which.max(fit.log$results$ROC)]
    # performance_training[7, 7] <- fit.log$results$F1[which.max(fit.log$results$ROC)]
    # if(length(unilabels) == 2){
    #   performance_training[8, 7] <- fit.log$results$Balanced_Accu[which.max(fit.log$results$ROC)]
    # }



    if(nrow(irisTest) > 0){

      logClasses <- predict(fit.log, newdata = irisTest, type = "prob")
      logClasses1 <- predict(fit.log, newdata = irisTest)
      logConfusion = confusionMatrix(data = logClasses1, irisTest$subtype)
      log.ROC <- roc(predictor = logClasses[,1], response = irisTest$subtype,
                     levels = rev(levels(irisTest$subtype)))
      performance_testing[1, 7] = as.numeric(log.ROC$auc)
      performance_testing[2, 7] = logConfusion$byClass[1]
      performance_testing[3, 7] = logConfusion$byClass[2]
      performance_testing[4, 7] = logConfusion$overall[1]
      performance_testing[5, 7] = logConfusion$byClass[5]
      performance_testing[6, 7] = logConfusion$byClass[6]
      performance_testing[7, 7] = logConfusion$byClass[7]
      if(length(unilabels) == 2){
        performance_testing[8, 7] <- logConfusion$byClass[11]
      }


    }


    #Deep learning#######

    # library(h2o)
    localH2O = h2o.init()
    prostate.hex <- as.h2o(irisTrain, destination_frame = "train.hex")
    hyper_params <- list(activation = c("Rectifier", "Tanh"),
                         hidden = list(c(100), c(200), c(10, 10), c(20, 20), c(50,50), c(30, 30, 30),
                                       c(25, 25, 25, 25)), input_dropout_ratio = c(0,0.05, 0.1),
                         l1 = seq(0, 1e-04, 1e-06), l2 = seq(0,1e-04, 1e-06),
                         train_samples_per_iteration = c(0,-2), epochs = c(50), momentum_start = c(0, 0.5),
                         rho = c(0.5, 0.99), quantile_alpha = c(0, 1), huber_alpha = seq(0,1))
    search_criteria = list(strategy = "RandomDiscrete", max_models = 100,
                           stopping_rounds = 5, stopping_tolerance = 0.01)


    dl_random_grid <- h2o.grid(algorithm = "deeplearning", grid_id = "dl_grid_randome1",
                               training_frame = prostate.hex, x = names(prostate.hex)[1:length(prostate.hex) -1],
                               y = "subtype", seed = 7, variable_importances = TRUE,
                               export_weights_and_biases = T, standardize = T, stopping_metric = "misclassification",
                               stopping_tolerance = 0.01, stopping_rounds = 2, score_duty_cycle = 0.025,
                               hyper_params = hyper_params, search_criteria = search_criteria,
                               nfolds = 10)


    grid <- h2o.getGrid("dl_grid_randome1", sort_by = "mse", decreasing = FALSE)
    grid@summary_table[1, ]
    best_model <- h2o.getModel(grid@model_ids[[1]])
    performance_training[1, 8] = as.numeric(best_model@model$cross_validation_metrics_summary['auc',1])
    performance_training[2, 8] = as.numeric(best_model@model$cross_validation_metrics_summary['recall',1])
    performance_training[3, 8] = as.numeric(best_model@model$cross_validation_metrics_summary['specificity',1])
    # performance_training[4, 8] = as.numeric(best_model@model$cross_validation_metrics_summary['accuracy',1])
    # performance_training[5, 8] = as.numeric(best_model@model$cross_validation_metrics_summary['precision',1])
    # performance_training[6, 8] = as.numeric(best_model@model$cross_validation_metrics_summary['recall',1])
    # performance_training[7, 8] = (2*performance_training[2, 8]*performance_training[5, 8])/(performance_training[2, 8]+performance_training[5, 8])
    #
    # if(length(unilabels) == 2){
    #   performance_training[8, 8]= (performance_training[2, 8]+performance_training[3, 8])/2}


      ##Get the result on testing data
      if(nrow(irisTest) > 0){

        perf = h2o.performance(best_model, as.h2o(irisTest, destination_frame = "test.hex"))
        performance_testing[1, 8] = as.numeric(h2o.auc(perf, h2o.find_threshold_by_max_metric(perf,
                                                                                              "f1"))[[1]])
        performance_testing[2, 8] = as.numeric(h2o.sensitivity(perf,
                                                               h2o.find_threshold_by_max_metric(perf, "f1"))[[1]])
        performance_testing[3, 8] = as.numeric(h2o.specificity(perf,
                                                               h2o.find_threshold_by_max_metric(perf, "f1"))[[1]])
        performance_testing[4, 8] = as.numeric(h2o.accuracy(perf,
                                                            h2o.find_threshold_by_max_metric(perf, "f1"))[[1]])
        performance_testing[5, 8] = as.numeric(h2o.precision(perf,
                                                             h2o.find_threshold_by_max_metric(perf, "f1")[[1]]))
        performance_testing[6, 8] = as.numeric(h2o.sensitivity(perf,
                                                               h2o.find_threshold_by_max_metric(perf, "f1"))[[1]])
        performance_testing[7, 8] = as.numeric(h2o.F1(perf, h2o.find_threshold_by_max_metric(perf,
                                                                                             "f1"))[[1]])
        if(length(unilabels) == 2){
          performance_testing[8, 8] = (performance_testing[2, 8] +
                                         performance_testing[3, 8])/2
        }
      }



      # performance_list[[k]] <- performance
      performance_testing_list[[k]] <- performance_testing
      performance_training_list[[k]] <- performance_training

      # performance_list <- list()
      if(length(unique(measurementLabels)) == 2){
        performance_training = matrix(rep(0, len = 24), nrow = 3)
        performance_testing = matrix(rep(0, len = 64), nrow = 8)
        # performance = matrix(rep(0, len = 64), nrow = 8)
      }else{
        performance_training = matrix(rep(0, len = 21), nrow = 3)
        performance_testing = matrix(rep(0, len = 64), nrow = 7)
        # performance = matrix(rep(0, len = 64), nrow = 7)
      }
    }

    list_test <- performance_testing_list
    list_train <- performance_training_list

    AUC_train <- lapply(list_train, function(x) x[1, ])
    AUC_test <- lapply(list_test, function(x) x[1, ])

    SENS_train <- lapply(list_train, function(x) x[2, ])
    SENS_test <- lapply(list_test, function(x) x[2, ])

    SPEC_train <- lapply(list_train, function(x) x[3, ])
    SPEC_test <- lapply(list_test, function(x) x[3, ])

    accuracy_test <- lapply(list_test, function(x) x[4, ])
    precision_test <- lapply(list_test, function(x) x[5, ])
    recall_test <- lapply(list_test, function(x) x[6, ])
    F1_test <- lapply(list_test, function(x) x[7, ])
    Balanced_accuracy_test <- lapply(list_test, function(x) x[8, ])

    output1 <- do.call(rbind, lapply(AUC_train, matrix, ncol = 8,
                                     byrow = TRUE))
    output2 <- do.call(rbind, lapply(AUC_test, matrix, ncol = 8,
                                     byrow = TRUE))
    output3 <- do.call(rbind, lapply(SENS_train, matrix, ncol = 8,
                                     byrow = TRUE))
    output4 <- do.call(rbind, lapply(SENS_test, matrix, ncol = 8,
                                     byrow = TRUE))
    output5 <- do.call(rbind, lapply(SPEC_train, matrix, ncol = 8,
                                     byrow = TRUE))
    output6 <- do.call(rbind, lapply(SPEC_test, matrix, ncol = 8,
                                     byrow = TRUE))

    output7 <- do.call(rbind, lapply(accuracy_test, matrix, ncol = 8,
                                     byrow = TRUE))
    output8 <- do.call(rbind, lapply(precision_test,
                                     matrix, ncol = 8, byrow = TRUE))
    output9 <- do.call(rbind, lapply(recall_test,
                                     matrix, ncol = 8, byrow = TRUE))
    output10 <- do.call(rbind, lapply(F1_test, matrix, ncol = 8,
                                     byrow = TRUE))
    output11 <- do.call(rbind, lapply(Balanced_accuracy_test,
                                     matrix, ncol = 8, byrow = TRUE))

    AUC_train_mean <- apply(output1, 2, mean)
    AUC_train_sd <- apply(output1, 2, sd)
    AUC_test_mean <- apply(output2, 2, mean)
    AUC_test_sd <- apply(output2, 2, sd)
    # AUC <- data.frame(AUC = t(cbind(t(AUC_train_mean), t(AUC_test_mean))))

    SENS_train_mean <- apply(output3, 2, mean)
    SENS_train_sd <- apply(output3, 2, sd)
    SENS_test_mean <- apply(output4, 2, mean)
    SENS_test_sd <- apply(output4, 2, sd)
    # SENS <- data.frame(t(SENS_train_mean), t(SENS_test_mean))

    SPEC_train_mean=apply(output5,2,mean)
    SPEC_train_sd=apply(output5,2,sd)
    SPEC_test_mean=apply(output6,2,mean)
    SPEC_test_sd=apply(output6,2,sd)
    # SPEC=data.frame(SPEC=t(cbind(t(SPEC_train_mean),t(SPEC_test_mean))))

    accuracy_test_mean <- apply(output7, 2, mean)
    accuracy_test_sd <- apply(output7, 2, sd)
    # accuracy <- data.frame(Accuracy=t(cbind(t(accuracy_test_mean))))

    precision_test_mean <- apply(output8, 2, mean)
    precision_test_sd <- apply(output8, 2, sd)
    # precision <- data.frame(Precision = t(cbind(t(precision_test_mean))))

    recall_test_mean <- apply(output9, 2, mean)
    recall_test_sd <- apply(output9, 2, sd)
    # recall <- data.frame(Recall=t(cbind(t(F1_test_mean))))

    F1_test_mean <- apply(output10, 2, mean)
    F1_test_sd <- apply(output10, 2, sd)
    # F1 <- data.frame(F1=t(cbind(t(F1_test_mean))))

    Balanced_accuracy_test_mean <- apply(output11, 2, mean)
    Balanced_accuracy_test_sd <- apply(output11, 2, sd)
    # Balanced_accuracy <- data.frame(Balanced_accuracy = t(cbind(t(Balanced_accuracy_test_mean))))

    trainingORtesting <- t(cbind(t(rep("training", 8)), t(rep("testing", 8))))
    testing_only <- t(t(rep("testing", 8)))
    training_only <- t(t(rep('training', 8)))

    performance_data_test <- data.frame(AUC = data.frame(AUC = t(c(t(AUC_test_mean), t(AUC_test_sd)))),
                                        SENS = data.frame(SENS = t((t(SENS_test_mean)))),
                                        SPEC = data.frame(SPEC = t((t(SPEC_test_mean)))),
                                        Accuracy = data.frame(Accuracy = t((t(accuracy_test_mean)))),
                                        Precision = data.frame(Precision = t((t(precision_test_mean)))),
                                        Recall = data.frame(Recall = t((t(recall_test_mean)))),
                                        F1 = data.frame(F1 = t((t(F1_test_mean)))),
                                        Balanced_accuracy = data.frame(Balanced_accuracy = t((t(Balanced_accuracy_test_mean)))),
                                        testing_only,
                                        Algorithm = (rep(t(c("RPART", "LDA", "SVM", "RF", "GBM", "PAM", "LOG", "DL")), 1)))

    performance_data_train <- data.frame(AUC = data.frame(AUC = t((t(AUC_train_mean)))),
                                         SENS = data.frame(SENS = t((t(SENS_train_mean)))),
                                         SPEC = data.frame(SPEC = t((t(SPEC_train_mean)))),
                                         training_only,
                                         Algorithm = (rep(t(c('RPART', 'LDA', 'SVM', 'RF', 'GBM', 'PAM', 'LoG', 'DL')), 1)))

    # library(reshape)
    melted_performance_data_test <- suppressMessages(melt(performance_data_test))
    melted_performance_data_train <- suppressMessages(melt(performance_data_train))

    sd_train <- data.frame(sd = rbind(t(t(AUC_train_sd)),t(t(SENS_train_sd)),t(t(SPEC_train_sd))))
    sd_test <- data.frame(sd = rbind(t(t(AUC_test_sd)),t(t(SENS_test_sd)),t(t(SPEC_test_sd)),
                                     t(t(accuracy_test_sd)),t(t(precision_test_sd)),t(t(recall_test_sd)),
                                     t(t(F1_test_sd)),t(t(Balanced_accuracy_test_sd))))
    melted_performance_data_train$sd <- sd_train
    melted_performance_data_test$sd <- sd_test

    # names(melted_performance_data_train)[1] <-
    #   names(melted_performance_data_test)[1] <- 'trainingORtesting'
    #
    # melted_performance_data <- rbind(melted_performance_data_train,
    #                                  melted_performance_data_test)

    # Algorithm <- NULL
    # value <- NULL
    # variable <- NULL
    # textLabels <- geom_text(aes(x = Algorithm, label = round(value, 2), fill = variable),
    #                         position = position_dodge(width = 1),
    #                         vjust = -0.5, size = 2)

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
    print(p1)

    if(nrow(irisTest) > 0){

      p2 <- ggplot(data = melted_performance_data_test,
                   aes(x = Algorithm, y = value, fill = variable)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                      position=position_dodge(.9)) +
        xlab("") + ylab("") + ggtitle("Testing") +
        theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 15, face = "bold"),
              axis.title = element_text(size = 14, face = "bold")) + labs(fill = "")
      print(p2)

    }




    plots <- list()
    performance <- list()

    performance$TrainingPerformance <- performance_data_train
    performance$TestingPerformance <- performance_data_test

    plots$TrainingPlot <- p1
    plots$TestingPlot <- p2


    return(c(performance, plots))
  }

