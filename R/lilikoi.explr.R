#' Exploratory analysis
#'
#' Performs source of variation test and build PCA and t-SNE plots to visualize important information.
#'
#' @param data is a input data frame for analysis with sample ids as row names and metabolite names or pathway names as column names.
#' @param demo.data is a demographic data frame with sample ids as row names, sample groups and demographic variable names as column names.
#' @param pca if TRUE, PCA plot will be out.
#' @param tsne if TRUE, T-SNE plot will be out.
#' @return Source of variation test results and PCA and t-SNE plot
#' @import parallel
#' @importFrom M3C pca tsne
#' @importFrom stats as.formula
#' @importFrom car Anova
#' @export
#' @examples
#' \donttest{
#' lilikoi.explr(data, demo.data, pca=TRUE, tsne=FALSE)
#' }


lilikoi.explr <- function(data, demo.data, pca=FALSE, tsne=FALSE){


  sovdata <- data
  pd <- demo.data

  # Calculate F
  sovdat <- sovdat[row.names(pd),]

  Ftab <- data.frame(Sample_Group = numeric(), stringsAsFactors = FALSE)
  for(i in 2:ncol(pd)){
    varname <- names(pd)[i]
    Ftab[varname] <- numeric()
  }


  calF <- function(probe = probecol){
    # probe = sovdata

    newdata <- pd
    pdnames <- names(newdata)
    newdata$beta <- probe

    formstr <- paste0(pdnames, collapse = ' + ')
    formstr <- paste0('beta ~ ', formstr)
    formstr <- as.formula(formstr)

    fit <- lm(formstr, data = newdata)

    aovfit <- Anova(fit, type = 3, singular.ok = TRUE)

    F <- aovfit$`F value`

    F <- F[2:(length(F)-1)]
    names(F) <- pdnames
    F <- as.data.frame(F, stringsAsFactors = FALSE)
    F <- as.data.frame(t(F))
    row.names(F) <- 1

    Ftab <- rbind(Ftab, F)

    return(Ftab)
  }


  sovdatlist <- list()
  for(i in 1:ncol(sovdat)){
    sovdatlist[[i]] <- sovdat[,i]
  }

  Ftab <- mclapply(X = sovdat, FUN = calF, mc.cores = 1)

  Ftab <- do.call(rbind, Ftab)

  Fmean <- colMeans(Ftab)

  Fmean <- Fmean[order(-Fmean)]

  Fmean <- data.frame(Factor = names(Fmean), Fstat = as.vector(Fmean), stringsAsFactors = FALSE)

  finalvars <- unique(c('Sample_Group', Fmean$Factor[Fmean$Fstat > 1]))


  sovplot <- function(restab = Fmean, clustername = 'Preeclampsia', plottype = 'MSS',
                      textsize = 20){

    resmean <- restab
    samplegroupidx <- match('Sample_Group', resmean$Factor)
    resmean$Factor[samplegroupidx] <- paste0(clustername, '_Control')

    if(plottype == 'MSS'){
      ytitle <- 'Mean Square'
      resmean <- resmean[order(-resmean$MSSstat),]
      resmean$Factor <- factor(resmean$Factor, levels = resmean$Factor, ordered = TRUE)

      p <- ggplot(data = resmean, mapping = aes(x = Factor, y = MSSstat, fill = Factor))
      print(
        p + geom_bar(stat = 'identity') +
          ggtitle('Source of Variance (Type 3 Anova)') +
          ylab(ytitle) +
          xlab('') +
          scale_fill_discrete(guide = FALSE) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 60, hjust = 1, size = textsize)) +
          theme(axis.text.y = element_text(size = textsize))
      )

    }else if(plottype == 'pval'){
      ytitle <- '-log2(p-val)'
      resmean <- resmean[order(-resmean$logpval),]
      resmean$Factor <- factor(resmean$Factor, levels = resmean$Factor, ordered = TRUE)

      p <- ggplot(data = resmean, mapping = aes(x = Factor, y = logpval, fill = Factor))
      print(
        p + geom_bar(stat = 'identity') +
          ggtitle('Source of Variance (Type 3 Anova)') +
          ylab(ytitle) +
          xlab('') +
          scale_fill_discrete(guide = FALSE) +
          geom_hline(yintercept = -log2(0.05), color = 'red', size = 1) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 60, hjust = 1, size = textsize)) +
          theme(axis.text.y = element_text(size = textsize))
      )
    }else{
      ytitle <- 'F statistic'
      resmean <- resmean[order(-resmean$Fstat),]
      resmean$Factor <- factor(resmean$Factor, levels = resmean$Factor, ordered = TRUE)

      p <- ggplot(data = resmean, mapping = aes(x = Factor, y = Fstat, fill = Factor))
      print(
        p + geom_bar(stat = 'identity') +
          ggtitle('Source of Variance (Type 3 Anova)') +
          ylab(ytitle) +
          xlab('') +
          scale_fill_discrete(guide = FALSE) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 60, hjust = 1, size = textsize)) +
          theme(axis.text.y = element_text(size = textsize))
      )
    }

  }

  plot = sovplot(restab = Fmean, plottype = 'F', textsize = 15)

  # return(c(Fmean,plot))

  # print(plot)

  returnList <- list()
  returnList$sovplot <- plot
  returnList$Fmean <- Fmean

  # return(Fmean)

  if (pca == TRUE){
    pca.plot <- pca(data)
    returnList$pca.plot <- pca.plot
  }

  if (tsne == TRUE){
    tsne.plot <- tsne(data)
    returnList$tsne.plot <- tsne.plot
  }

  return(returnList)

}







