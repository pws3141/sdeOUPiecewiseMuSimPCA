#! /usr/bin/env Rscript

# Simulate loadings using OU process
# with piecewise cont mu
# dX = gamma (mu(t) - X) dt + dW

#simTime <- as.integer(args)
simTime <- 100000

cat("Initialising...\n")
set.seed(10)

library(MASS) # for mvrnorm
suppressMessages(library(ggplot2))
library(Cairo) # for cairo type in ggsave
library(rlist)

source("../sdeParamEstimationOUBishwal/sdeEMSchemeOU.R")

cat("loading PCA maximum likelihood estimators...\n")
# load loadings to get initial points
mslpPCALoadings <- as.matrix(read.csv("../../meetings/19.06.17/hadcmPCA/loadings10.csv", 
                                      header = FALSE))
mslpPCALoadingsInitial <- mslpPCALoadings[1,]

mleEstimates <- list.load(file = "../../meetings/19.08.01/pcaOUPiecewiseMu/mleEstimates.yaml")

cat("simulating one dimensional SDE from MLE estimate...\n")
# now simulate points for each month using muMLE and varMLE
numPointsMonth <- floor(simTime / 12)

# define t: not sure what delta t to use

t <- seq(from = 0, by = 0.001, length = simTime)
#t <- seq(from = 0, to = simTime - 1, by = 1)

mle <- mleEstimates[[1]]
x0 <- mslpPCALoadingsInitial[1]
mut <- mle$muHat
gammat <- mle$gammaHat
path <- ouProcessPiecewiseEMApproximation(x0 = x0, 
                         gammat = gammat, mut = mut,
                         sigma = 1, tau = t) 

cat("plotting sde..\n")
# plot the first simulated loadings
# matrix of one simulation for each loading

whichMu <- (floor(t) %% 12) + 1
muTPlot <- mut[whichMu] 

X <- list(mu = muTPlot, t = path$t, y = path$y)
dfX <- as.data.frame(X)
names(dfX) <- c("mu", "time", "X")

plotName <- "simulatedSDE.pdf"
gg <- ggplot(dfX, aes(x = time, y = X)) + 
        geom_path() +
        geom_point(color = "cyan", size = 0.4) +
        geom_path(aes(x = time, y = mu), 
                  color = "blue") +
        labs(x = "time",
             y = "Y") #+
        #coord_fixed()
#
ggsave(plotName, plot = gg, 
       width = 30, height = 15, units = "cm")

cat("done...\n")



