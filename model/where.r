#!/usr/bin/Rscript
library(rEDM)
library(zoo)
source( "../helpers/helper.r" )
source( "../helpers/plotting.r" )

args <- commandArgs(trailingOnly = TRUE)

if( !is.na(args[1]) ) {
    print("You sure you wanna erase everything and re-download?" )
    answer <- readLines( "stdin", n=1 )
    
    if( answer == "yes" ) {
        print( "OK. Please insert password.")
        system("rm -f ~/edm/weight-pred/model/run/*")
        system( "scp -r daon@access.cims.nyu.edu:~/weight-pred/model/runs/ ." )
    }
}

load( "runs/parameters.Rdata" )

## Get the true values we were trying to predict
df <- read.csv("originals/three_species.csv",
                   header = TRUE,
                   sep = "," )
df <- lag_every_variable(df, n_lags)
truth <- df$x_p1[pred[1]:pred[2]]

weighted_pred <- read.table("runs/weighted_predictions_x.csv",
                            na = "NA" )
                            ## row.names = FALSE,
                            ## col.names = FALSE)

mve_pred <- read.table("runs/mve_predictions_x.csv",
                       na = "NA" )


weighted_rho <- cor( t(weighted_pred), truth, use = "complete.obs" )
weighted_rho <- colMeans( matrix(weighted_rho, nrow = n_samp) )

mve_rho <- cor( t(mve_pred), truth, use = "complete.obs" )
mve_rho <- colMeans( matrix(mve_rho, nrow = n_samp) )

x11()
plot(lib_sizes,
     weighted_rho,
     type = "l",
     col = "black" )
lines(lib_sizes,
      mve_rho,
      col = "red" )
hold(1)

## libs <- read.table("runs/libraries_x.txt",
##                    header = FALSE,
##                    sep = "\t" )
