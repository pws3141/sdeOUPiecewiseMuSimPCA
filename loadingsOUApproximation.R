#! /usr/bin/env Rscript

# Assume data simulated from OU process
# dX = gamma (mu(t) - X) dt + dW
# with mu(t) piecewise constant function, mod(12)

cat("Initialising...\n")

#library(jvcoords)
#library(chron)
#library(RColorBrewer)
#library(lattice)
#library(ncdf4)
#library(reshape2) #for melt
suppressMessages(library(ggplot2))
library(rlist)
#library(rgdal) # reading shapefile data
#library(scales) # for rescale in ggplot
#library(Cairo) # for cairo type in ggsave

source("../sdeParamEstimationOUBishwal/mleOU.R")

cat("loading PCA data...\n")
mslpPCALoadings <- as.matrix(read.csv("../../19.06.17/hadcmPCA/loadings10.csv", 
                                      header = FALSE))
numLoadings <- ncol(mslpPCALoadings)
timeSteps <- nrow(mslpPCALoadings)


dfX <- as.data.frame(mslpPCALoadings)
names(dfX) <- paste0("X", seq_len(numLoadings))

# define t: not sure what delta t to use
t <- seq(from = 0, to = timeSteps - 1, by = 1)
splitMonth <- lapply(seq_len(12), function(i) {
                             seq(from = i, to = timeSteps, by = 12)
                })

cat("finding MLE...\n")

# initialise matrix of mle estimates
mleEstimates <- vector(mode = "list", length = numLoadings)
names(mleEstimates) <- paste0("loading", 1:numLoadings)
#mleEstimates <- lapply(mleEstimates, function(x) {
                               #mat <- matrix(ncol = 2, nrow = 12)
                               #rownames(mat) <- month.abb
                               #colnames(mat) <- c("gammaHat", "muHat")
                               #mat
                #})

for(i in 1:numLoadings) {
        XTmp <- dfX[[i]]
        XTmp <- XTmp[1:timeSteps]
        # split time series in 12
        XSplit <- lapply(splitMonth, function(i) XTmp[i])
        mleMuTmp <- sapply(XSplit, function(x) {
                                 tTmp <- 1:length(x) - 1
                                 mleOU(X = x, t = tTmp)
                              })
        colnames(mleMuTmp) <- month.abb
        mleGammaTmp <- mleOU(X = XTmp, t = t)
        mleEstimates[[i]] <- list(muHat = t(mleMuTmp[2, ,drop=FALSE]),
                                  gammaHat = mleGammaTmp$gammaHat)
        #mleEstimates[[i]]$muHat <- t(mleMuTmp[2, ,drop=FALSE])
        #mleEstimates[[i]]$gammaHat <- mleGammaTmp$gammaHat
}

fname <- "mleEstimates.yaml"
cat("writing ", fname, "...\n")
list.save(mleEstimates, fname)

cat("done...\n")

