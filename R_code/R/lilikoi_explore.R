#' Exploratory analysis
#'
#' Perform source of variation test and build t-SNE plot.
#'
#' @param sovdata input data for analysis
#' @param func regression function
#' @param selected_Pathways_Weka Selected top pathways from the featureSelection function
#' @param Metabolite_pathway_table Metabolites mapping table
#' @return Source of analysis test results and t-SNE plot
#' @export
#' @examples lilikoi_explore(sovdata = sovdata, func = function)


#Calculate F#####
# Ftab <- data.frame(Sample_Group = numeric(), stringsAsFactors = FALSE)
# for(i in 2:ncol(pd)){
#   varname <- names(pd)[i]
#   Ftab[varname] <- numeric()
# }
#
# # probe = sovdat
# calF <- function(probe = probecol){
#   library(car)
#
#   newdata <- pd
#   pdnames <- names(newdata)
#   newdata$beta <- probe
#
#   formstr <- paste0(pdnames, collapse = ' + ')
#   formstr <- paste0('beta ~ ', formstr)
#   formstr <- as.formula(formstr)
#
#   fit <- lm(formstr, data = newdata)
#
#   aovfit <- Anova(fit, type = 3, singular.ok = TRUE)
#
#
#   F <- aovfit$`F value`
#
#   F <- F[2:(length(F)-1)]
#   names(F) <- pdnames
#   F <- as.data.frame(F, stringsAsFactors = FALSE)
#   F <- as.data.frame(t(F))
#   row.names(F) <- 1
#
#
#   Ftab <- rbind(Ftab, F)
#
#   return(Ftab)
# }
#
#
# library(parallel)
#
# sovdatlist <- list()
#
# for(i in 1:ncol(sovdat)){
#   sovdatlist[[i]] <- sovdat[,i]
# }
#
#
# Ftab <- mclapply(X = sovdat, FUN = calF, mc.cores = 1)
# #Ftab <- apply(X = sovdat[,1:50], MARGIN = 2, FUN = calF)
#
# Ftab <- do.call(rbind, Ftab)
#
#
# Fmean <- colMeans(Ftab)
#
# Fmean <- Fmean[order(-Fmean)]
#
# Fmean <- data.frame(Factor = names(Fmean), Fstat = as.vector(Fmean), stringsAsFactors = FALSE)
#
# finalvars <- unique(c('Sample_Group', Fmean$Factor[Fmean$Fstat > 1]))
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#






