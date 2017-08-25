#!/usr/bin/Rscript
library(rEDM)
source( "helper.r" )
source( "ccm.r" )

df <- read.csv( "data/univariate_data.csv", header = TRUE )
simplex_output <- simplex( df$chl, E =1:15 )

x11()
plot(simplex_output$E, simplex_output$rho, type = "l", xlab = "Embedding Dimension (E)", 
    ylab = "Forecast Skill (rho)")
hold(0)


show_ccm <- function( var, tp ) {
    print( paste0( "Chl X ", var, " tp = ", tp, ", rho = ", chl_xmap_var( var, tp = -2 )$rho ) )
    print( paste0(  var, " X Chl, tp = ", tp, ", rho = ", var_xmap_chl( var, tp = -2 )$rho ) )
}

show_ccm( "nitrate" , -2 )
show_ccm( "silicate", -2 )


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
