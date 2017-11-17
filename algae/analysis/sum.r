#!/usr/bin/Rscript
library(rEDM)
source( "helpers/helper.r" )
source( "helpers/plotting.r" )





## Load the data: df, chl_threshold and norm_threshold.
## args <- commandArgs(trailingOnly = TRUE)
## if( length(args) == 0 ) {
##     xtension <- "full"
## } else if( grepl("bloom", args, ignore.case = TRUE ) ) {
##     xtension <- "bloom"
## } else {
##     stop( "Choose bloom or full!" )
## }

load( "data/leftover_predictions_full_loocv.Rdata" )
load( "data/processed_block.Rdata" )

obs_cor <- function(x) return( cor (    x+main_pred , df$chl_p1wk,  use = "complete.obs" ) )
obs_mae <- function(x) return( mean(abs(x+main_pred - df$chl_p1wk), na.rm = TRUE         ) )

rhos <- apply(leftover_pred, 1, obs_cor )
maes <- apply(leftover_pred, 1, obs_mae )
score <- rank(rhos) + rank(-maes)

ind  <- which.max( score )
best <- paste(other_vars[ combinations[ , ind ] ], collapse = ", " )
print( paste0( "Best model (comb) with rho = ", rhos[ind], " and MAE = ", maes[ind] ) )
print( best )

ind  <- which.max( rhos )
best <- paste(other_vars[ combinations[ , ind ] ], collapse = ", " ) 
print( paste0( "Best model (rho)  with rho = ", rhos[ind], " and MAE = ", maes[ind] ) )
print( best )

ind <- which.max( -maes )
best <- paste(other_vars[ combinations[ , ind ] ], collapse = ", " ) 
print( paste0( "Best model (MAE)  with rho = ", rhos[ind], " and MAE = ", maes[ind] ) )
print( best )

model1 <- c( "chl", other_vars[ combinations[ind] ] )
model2 <- c( "chl", "silicate", "silicate_1wk", "AvgDens" )

