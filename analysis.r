#!/usr/bin/Rscript
library(rEDM)
source( "plotting.r" )
source( "helper.r" )
source( "ccm.r" )

## df <- read.csv( "data/univariate_data.csv", header = TRUE )
## simplex_output <- simplex( df$chl, E =1:15 )

## x11()
## plot(simplex_output$E,
##      simplex_output$rho,
##      type = "l",
##      xlab = "Embedding Dimension (E)", 
##      ylab = "Forecast Skill (rho)" )
## hold(0)

libsizes  <- (1:8) * 10
n_samples <- 20
show_ccm( "nitrate" , -2, libsizes, n_samples, 4 )
## show_ccm( "silicate", -2, libsizes, n_samples, 4 )


## x11()
## plot(E,
##      rhos,
##      type = "l",
##      col = "blue",
##      main = paste0("Prediction of ", predictee, " using ", predictor, " data using ", numLags, " lags") )
## hold(1)

## perm <- order(rho)
## print( paste0("Lags = ",
##               E[length(E)],
##               ", max rho = ",
##               rho[perm[length(rho)]],
##               ", E = ",
##               perm[length(rho)],
##               ", 2nd rho = ",
##               rho[perm[length(rho)-1]],
##               ", E = ",
##               perm[length(rho)-1] )
##       )
