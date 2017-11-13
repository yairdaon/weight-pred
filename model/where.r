#!/usr/bin/Rscript
library(rEDM)
library(zoo)
source( "../helpers/helper.r" )
source( "../helpers/plotting.r" )

results <- function( method, variable )
    {
        load( paste0("runs/", method, "_parameters.Rdata" ) )
        
        ## Get the true values we were trying to predict
        df <- read.csv(filename,
                       header = TRUE,
                       sep = "," )
        df <- lag_every_variable(df, n_lags)
        truth <- df[pred[1]:pred[2], paste0( variable, "_p1" ) ]

        predictions <- read.table(paste0("runs/", method, "_predictions_", variable, ".csv"),
                                  na = "NA",
                                  sep =",") 

        ## Calculate mean prediction skill per library size
        rhos <- mean_cor( predictions, truth, nrow = n_samp)
        return( rhos ) 

}

args <- commandArgs( trailingOnly = TRUE )
if( length(args) > 0 ) {
    
    print("Erase everything and re-download? Press Enter to skip, yes to download." )
    answer <- readLines( "stdin", n=1 )
    
    if( answer == "yes" ) {
        print( "OK. Please insert password.")
        system("rm -f ~/edm/weight-pred/model/run/*")
        system( "scp -r daon@access.cims.nyu.edu:~/weight-pred/model/runs/ ." )
    }
}

rhos <- results( "ml", "y" )
print( "Max Likelihood" )
print( rhos )

rhos <- results( "exp", "y" )
print( "Exponential Weights" )
print( rhos )

rhos <- results( "mve", "y" )
print( "MVE" )
print( rhos )

## x11()
## plot(lib_sizes,
##      weighted_rho,
##      type = "l",
##      col = "black" )
## lines(lib_sizes,
##       mve_rho,
##       col = "red" )
## hold(1)

## libs <- read.table("runs/libraries_x.txt",
##                    header = FALSE,
##                    sep = "\t" )
