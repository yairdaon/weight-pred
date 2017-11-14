#!/usr/bin/Rscript
library(rEDM)
library(zoo)
source( "helpers/helper.r" )
source( "helpers/plotting.r" )

results <- function( method, variable )
{
    load( paste0("model/runs/", method, "_parameters.Rdata" ) )
    
    ## Get the true values we were trying to predict
    df <- read.csv(filename,
                   header = TRUE,
                   sep = "," )
    df <- lag_every_variable(df, n_lags)
    truth <- df[pred[1]:pred[2], paste0( variable, "_p1" ) ]
    
    predictions <- read.table(paste0("model/runs/", method, "_predictions_", variable, ".csv"),
                              na = "NA",
                              sep =",") 
    
    ## Calculate mean prediction skill per library size
    rhos <- mean_cor( predictions, truth, nrow = n_samp)
    return( list( rhos = rhos, lib_sizes = lib_sizes ) )
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

variable <- "y"
minvar      <- results( "minvar",  variable )
exponential <- results( "exp",     variable )
uniform     <- results( "uniform", variable )

x11()
plot(minvar$lib_sizes,
     minvar$rhos,
     type = "l",
     col = "black",
     xlab = "Library size",
     ylab = "Prediction skill",
     ylim = c(0.5,1)
     )
lines(exponential$lib_sizes,
      exponential$rhos,
      col = "red" )
lines(uniform$lib_sizes,
      uniform$rhos,
      col = "blue" )
title(paste0("Predicting ", variable) )
legend( x = "topleft",
       legend = c( "Min-Var weights", "Exp(-Var)", "Uniform MVE" ),
       col = c( "black", "red", "blue" ),
       lwd = 1)
hold(1)

## libs <- read.table("runs/libraries_x.txt",
##                    header = FALSE,
##                    sep = "\t" )
