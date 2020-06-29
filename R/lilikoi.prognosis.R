#' Pathway-based prognosis model
#'
#' Fits a Cox proportional hazards regression model or a Cox neural network model
#' to predict survival results.
#'
#' @param event survival event
#' @param time survival time
#' @param exprdata dataset for penalization, with id in the rownames and pathway or metabolites
#' names in the column names.
#' @param percent train-test separation percentage
#' @param alpha denote which penalization method to use.
#' @param nfold fold number for cross validation
#' @param method determine the prognosis index, "quantile", "quantile" or "ratio".
#' @param cvlambda determine the lambda for prediction, "lambda.min" or "lambda.1se".
#' @param python.path saved path for python3
#' @param path saved path for the inst file in lilikoi2
#' @param coxnnet if TRUE, coxnnet will be used.
#' @param coxnnet_method the algorithm for gradient descent. Includes standard gradient descent ("gradient"), Nesterov accelerated gradient "nesterov" and momentum gradient descent ("momentum").
#' @import reticulate survminer
#' @importFrom stats complete.cases median predict quantile
#' @importFrom utils write.table
#' @importFrom glmnet cv.glmnet
#' @importFrom survival Surv survfit survdiff coxph
#' @return A list of components: \item{c_index}{C-index of the Cox-PH model} \item{difftest}{Test results of the survival curve difference test}
#' \item{survp}{Kaplan Meier plot}
#' @export
#' @examples
#' \donttest{
#' inst.path = path.package('lilikoi2', quiet = FALSE) # path = "lilikoi/inst/", use R to run
#' lilikoi.prognosis(event, time, exprdata, percent=NULL, alpha=0, nfold=5, method="median",
#'   cvlambda=NULL,python.path=NULL, inst.path=inst.path,coxnnet=FALSE,coxnnet_method="gradient")
#' }


lilikoi.prognosis <- function(event, time, exprdata, percent=NULL, alpha=0,
                      nfold=5, method="median", cvlambda=NULL,
                      python.path=NULL, path=NULL,coxnnet=FALSE,
                      coxnnet_method="gradient"){

  if (coxnnet==FALSE){
    if (is.null(percent)){
      s <- Surv(time=time,event=event)

      cvfit <- cv.glmnet(exprdata, s, family = "cox", maxit = 1000, alpha=alpha, nfolds = nfold)
      score <- predict(cvfit, newx = exprdata, s=cvlambda)

      PIscore <- as.data.frame(score, drop = F)
      colnames(PIscore)  <- "score"
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

      s <- Surv(time=time_train,event=event_train)

      cvfit <- cv.glmnet(exprdata_train, s, family = "cox", maxit = 1000, alpha=alpha)
      score <- predict(cvfit, newx = exprdata_test, s=cvlambda)

      PIscore <- as.data.frame(score, drop = F)
      colnames(PIscore)  <- "score"
      PIscore$EVENT <- event_test
      PIscore$TIME <- time_test
    }
  }

  if(coxnnet==TRUE){

    if (is.null(percent)){
      original.path <- getwd() # Get original path
      on.exit(setwd(original.path))

      setwd(path)
      on.exit(setwd(original.path))

      write.table(exprdata, 'x.csv', row.names=FALSE, col.names=FALSE, sep = ",")
      write.table(time, 'ytime.csv', row.names=FALSE, col.names=FALSE, sep = ",")
      write.table(event, 'ystatus.csv', row.names=FALSE, col.names=FALSE, sep = ",")


      use_python(python.path) # user input: terminal which python3
      # repl_python() # Install miniconda; then exit.
      # exit
      use_virtualenv("r-reticulate")
      py_available(TRUE)

      # py_install("sklearn",pip = TRUE) # Skip the installation if already done so.
      # py_install("numpy",pip = TRUE)
      # py_install("theano",pip = TRUE)

      # Read in python code
      # pathnew = file.path(path, 'cox_nnet3')
      source_python("L2cross_nopercent.py")
      result <- L2cross_nopercent(path,coxnnet_method,percent,nfold)

      # Read in the results on the test set
      as.numeric(readLines("theta.csv")) -> score
      as.numeric(readLines("ytime.csv")) -> time
      as.numeric(readLines("ystatus.csv")) -> status


      setwd(original.path)
      on.exit(setwd(original.path))

      PIscore <- as.data.frame(score, drop = F)
      colnames(PIscore)  <- "score"
      PIscore$EVENT <- status
      PIscore$TIME <- time
    }else{
      original.path <- getwd() # Get original path
      on.exit(setwd(original.path))

      setwd(path)
      on.exit(setwd(original.path))

      write.table(exprdata, "x.csv", row.names=FALSE, col.names=FALSE, sep = ",")
      write.table(time, "ytime.csv", row.names=FALSE, col.names=FALSE, sep = ",")
      write.table(event, "ystatus.csv", row.names=FALSE, col.names=FALSE, sep = ",")


      use_python(python.path) # user input: terminal which python3
      # repl_python() # Install miniconda; then exit.
      # exit
      use_virtualenv("r-reticulate")
      py_available(TRUE)

      # py_install("sklearn",pip = TRUE) # Skip the installation if already done so.
      # py_install("numpy",pip = TRUE)
      # py_install("theano",pip = TRUE)

      # Read in python code
      # pathnew = file.path(path, 'cox_nnet3')
      source_python("L2cross.py")
      result <- L2cross(path,coxnnet_method,percent,nfold)

      # Read in the results on the test set
      as.numeric(readLines("theta.csv")) -> score
      as.numeric(readLines("ytime_test.csv")) -> time
      as.numeric(readLines("ystatus_test.csv")) -> status

      setwd(original.path)
      on.exit(setwd(original.path))

      PIscore <- as.data.frame(score, drop = F)
      colnames(PIscore)  <- "score"
      PIscore$EVENT <- status
      PIscore$TIME <- time
    }

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

  if (method == "ratio"){
    PI = log(sum(event)/(length(event)-sum(event)))
    for(i in 1:nrow(PIscore))
    {
      if(PIscore[i,"score"]<=PI)
        PIscore[i,"highlabel"]<-0

      else if (PIscore[i,"score"]>PI)
        PIscore[i,"highlabel"]<-1
    }
  }

  PIscore_na <- PIscore[complete.cases(PIscore),]

  # attach(PIscore_na)
  surv <- Surv(PIscore_na$TIME, PIscore_na$EVENT)
  sum.surv <- summary(coxph(surv ~ PIscore_na$highlabel))
  c_index <- sum.surv$concordance

  # print(c_index)

  difftest <- survdiff(Surv(TIME, EVENT) ~ highlabel, data=PIscore_na)

  # print(difftest)

  fit_surv = surv_fit(Surv(PIscore_na$TIME, EVENT) ~ highlabel, data = PIscore_na)
  survp <- ggsurvplot(fit_surv, data=PIscore_na)

  # return(survp)

  returnList <- list("c_index"=c_index, "difftest"=difftest, "survp"=survp)
  return(returnList)

}





