#! /usr/bin/env Rscript

# Simulate loadings using OU process
# with piecewise cont mu
# dX = gamma (mu(t) - X) dt + dW

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
numSims <- as.integer(args)
#simTime <- as.integer(args)
simTime <- 6000

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

mleEstimates <- list.load(file = "mleEstimates.yaml")
numLoadings <- length(mleEstimates)

cat("simulating loadings from MLE estimates...\n")
# now simulate points for each month using muMLE and varMLE
numPointsMonth <- floor(simTime / 12)

# define t: not sure what delta t to use

t <- seq(from = 0, by = 0.01, length = simTime)
tSwitch <- 1 # switch param in EM approx function
#t <- seq(from = 0, to = simTime - 1, by = 1)

simulatedLoadings <- vector(mode = "list", length = numLoadings)
simulatedLoadings <- lapply(seq_len(numLoadings), 
                        function(i) {
                                tmp <- matrix(ncol = numSims, nrow = simTime)
                                colnames(tmp) <- paste0("simulation", 1:numSims)
                                tmp
                        })
names(simulatedLoadings) <- paste0("loading", 1:numLoadings)

simulatedAverageLoadings <- vector(mode = "list", length = numLoadings)
simulatedAverageLoadings <- lapply(seq_len(numLoadings), 
                        function(i) {
                                tmp <- matrix(ncol = numSims, nrow = floor(max(t) + 1))
                                colnames(tmp) <- paste0("simulation", 1:numSims)
                                tmp
                        })
names(simulatedAverageLoadings) <- paste0("loading", 1:numLoadings)

for (j in 1:numSims) {
        simulatedLoadingsTmp <- lapply(seq_len(numLoadings), function(i) {
               mle <- mleEstimates[[i]]
               x0 <- mslpPCALoadingsInitial[i]
               mut <- mle$muHat
               gammat <- mle$gammaHat
               path <- ouProcessPiecewiseEMApproximation(x0 = x0, 
                                         gammat = gammat, mut = mut,
                                         sigma = 1, tau = t) 
               path$y
              })
        yMonthlyAverage <- lapply(simulatedLoadingsTmp, function(y) {
                              tmp <- sapply(seq_len(floor(max(t)) + 1), function(i) {
                                      tMonth <- which(t < i & t >= i - 1)
                                      yMonth <- y[tMonth]
                                      mean(yMonth)
                              })
                              tmp
              })

        simulatedLoadings <- lapply(seq_len(numLoadings), function(i) {
                                simTmp <- simulatedLoadingsTmp[[i]]
                                simulatedLoadings[[i]][, j] <- simTmp
                                simulatedLoadings[[i]]
                              })
        simulatedAverageLoadings <- lapply(seq_len(numLoadings), function(i) {
                                simTmp <- yMonthlyAverage[[i]]
                                simulatedAverageLoadings[[i]][, j] <- simTmp
                                simulatedAverageLoadings[[i]]
                              })
        
}

cat("plotting loadings...\n")
# plot the first simulated loadings
# matrix of one simulation for each loading
plotSimulatedLoadings <- lapply(simulatedLoadings, function(x) x[, 1])
plotSimulatedLoadings <- do.call(cbind, plotSimulatedLoadings)

dfX <- as.data.frame(plotSimulatedLoadings)
names(dfX) <- paste0("X", seq_len(numLoadings))

# plot the first simulated loadings
# matrix of one simulation for each loading

#

for (i in 1:4) {
        plotName <- paste0("PCASimulatedLoading", i, ".png")
        mle <- mleEstimates[[i]]
        muTmp <- mle$muHat
        whichMu <- (floor(t / tSwitch) %% length(muTmp)) + 1
        muTPlot <- muTmp[whichMu] 
        #muTPlot <- rep(mu, times = length(t))

        X <- list(mu = muTPlot, t = t, x = plotSimulatedLoadings[, i])
        dfXTmp <- as.data.frame(X)
        names(dfXTmp) <- c("mu", "time", "X")
        mainLabs <- paste0("Simulated Loading ", i, ": Gamma = ", mle$gammaHat)

        gg <- ggplot(dfXTmp, aes(x = time, y = X)) + 
                geom_path() +
                geom_point(color = "cyan", size = 0.4) +
                geom_path(aes(x = time, y = mu), 
                          color = "blue") +
                labs(title = mainLabs, x = "time",
                     y = "Y") #+
                #coord_fixed()
        #
        ggsave(plotName, plot = gg, path = "./simulatedLoadings/",
               width = 30, height = 15, units = "cm")

        for (j in (i+1):4) {
                if (i == 4) break
                mainLabs <- paste0("Simulated Loadings ", i, " and ", j, 
                                   ": Gamma = ", mle$gammaHat)
                plotName <- paste0("PCAtwoSimulatedLoadings", i, j, ".png")
                ii <- paste0("X", i)
                jj <- paste0("X", j)
                gg <- ggplot(dfX, aes_string(x = ii, y = jj)) + 
                        geom_path() +
                        geom_point(color = "cyan", size = 0.4) +
                        labs(title = mainLabs, 
                             x = paste0("Sim Loading", i), 
                             y = paste0("Sim Loading", j)) +
                        coord_fixed()
                #
                ggsave(plotName, plot = gg, path = "./simulatedLoadings/", 
                       width = 15, height = 15, units = "cm")
        }
}

plotSimulatedAverageLoadings <- lapply(simulatedAverageLoadings, function(x) x[, 1])
plotSimulatedAverageLoadings <- do.call(cbind, plotSimulatedAverageLoadings)

dfX <- as.data.frame(plotSimulatedAverageLoadings)
names(dfX) <- paste0("X", seq_len(numLoadings))

for (i in 1:3) {
        for (j in (i+1):4) {
                plotName <- paste0("PCAtwoSimulatedAverageLoadings", i, j, ".png")
                ii <- paste0("X", i)
                jj <- paste0("X", j)
                gg <- ggplot(dfX, aes_string(x = ii, y = jj)) + 
                        geom_path() +
                        geom_point(color = "cyan", size = 0.4) +
                        labs(x = paste0("Sim Loading", i), 
                             y = paste0("Sim Loading", j)) +
                        coord_fixed()
                #
                ggsave(plotName, plot = gg, path = "./simulatedLoadings/", 
                       width = 15, height = 15, units = "cm")
        }
}


fname <- "simulatedAverageLoadings.yaml"
cat("writing", fname, "...\n")
list.save(simulatedAverageLoadings, fname)


fname <- "simulatedLoadings.yaml"
cat("writing", fname, "...\n")
list.save(simulatedLoadings, fname)

cat("done...\n")


