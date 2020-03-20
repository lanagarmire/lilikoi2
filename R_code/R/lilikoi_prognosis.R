#' Pathway-based prognosis model
#'
#'
#' @param PDSmatrix Pathway deregulation score matrix`1àQDF`
#' @param event survival event
#' @param time survival time
#' @param exprdata dataset for penalization
#' @param percent train-test separation percentage
#' @param alpha denote which penalization method to use.
#' @param nfold fold number for cross validation
#' @param method determine the prognosis index
#' @return
#' @export
#' @examples lilikoi_prognosis(PDSmatrix, event, time, exprdata, percent, alpha, nfold, method)


# PDSmatrix,  event, time, alpha, percent, exprdata, nfold, method

# Input PDSmatrix,

lilikoi_prognosis <- function(PDSmatrix, event, time, exprdata, percent, alpha, nfold, method){
  if (is.null(percent)){

    library(survival)
    s <- Surv(time=time,event=event)

    library(penalized)
    if (alpha == 0){
      optjc <- optL2(s, penalized=exprdata, fold=nfold) # Ridge
      fit <- penalized(s, penalized=exprdata, lambda2=optjc$lambda)
    }

    if (alpha == 1){
      optjc <- optL1(s, penalized=exprdata, fold=nfold)
      fit <- penalized(s, penalized=exprdata, lambda1=optjc$lambda)
    }

    if (alpha == 0.5){
      # optjc <- optL1(s, penalized=exprdata, fold=nfold)
      fit <- penalized(s, penalized=exprdata, lambda1=1, lambda2=1)
    }

    # score <- predict(fit, exprdata) # predict
    score <- fitted(fit, exprdata)
    score <- log(score)
    PIscore <- as.data.frame(score)
    PIscore$ID <- rownames(exprdata_tumor)

    PIscore$EVENT <- event
    PIscore$TIME <- time

  }else{
    # Train/test split -- percent
    smp_size <- floor(percent * length(time))
    train_ind <- sample(seq_len(length(time)), size = smp_size)

    time_train <- time[train_ind]
    time_test <- time[-train_ind]

    event_train <- event[train_ind]
    event_test <- event[-train_ind]

    exprdata_train <- exprdata[train_ind,]
    exprdata_test <- exprdata[-train_ind,]

    library(survival)
    s <- Surv(time=time_train,event=event_train)

    library(penalized)
    if (alpha == 0){
      optjc <- optL2(s, penalized=exprdata_train, fold=nfold)
      fit <- penalized(s, penalized=exprdata_train, lambda2=optjc$lambda)
    }

    if (alpha == 1){
      optjc <- optL1(s, penalized=exprdata_train, fold=nfold)
      fit <- penalized(s, penalized=exprdata_train, lambda1=optjc$lambda)
    }

    if (alpha == 0.5){
      # optjc <- optL1(s, penalized=exprdata_train, fold=nfold)
      fit <- penalized(s, penalized=exprdata_train, lambda1=1, lambda2=1)
    }

    score <- predict(fit, exprdata_test) # predict
    score <- log(score)
    PIscore <- as.data.frame(score)
    PIscore$ID <- rownames(exprdata_tumor)

    PIscore$EVENT <- event_test
    PIscore$TIME <- time_test
  }


  # Predict high label
  if (method == "median"){
    PI = median(PIscore$score)
    for(i in 1:nrow(PIscore))
    {
      if(PIscore[i,"score"]<=PI)
        PIscore[i,"highlabel"]<-0

      else if (PIscore[i,"score"]>PI)
        PIscore[i,"highlabel"]<-1
    }
  }

  if (method == "quantile"){
    PI = quantile(PIscore$score)
    for(i in 1:nrow(PIscore))
    {
      if(PIscore[i,"score"]<=PI[[2]])
        PIscore[i,"highlabel"]<-0

      else if (PIscore[i,"score"]>PI[[4]])
        PIscore[i,"highlabel"]<-1
    }
  }

  PIscore_na <- PIscore[complete.cases(PIscore),]

  attach(PIscore_na)
  surv <- Surv(TIME, EVENT)
  sum.surv <- summary(coxph(surv ~ highlabel))
  c_index <- sum.surv$concordance

  print(c_index)

  difftest <- survdiff(Surv(TIME, EVENT) ~ highlabel)

  print(difftest)

  library(survminer)
  fit_surv = survfit(Surv(TIME, EVENT) ~ highlabel, data = PIscore_na)
  survp <- ggsurvplot(fit_surv, data = PIscore_na)

  print(survp)

}

