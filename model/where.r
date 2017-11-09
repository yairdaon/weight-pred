#!/usr/bin/Rscript
library(rEDM)
library(zoo)
source( "../helpers/helper.r" )
source( "../helpers/plotting.r" )

print("Erase everything and re-download? Press Enter to skip, yes to download." )
answer <- readLines( "stdin", n=1 )

if( answer == "yes" ) {
    print( "OK. Please insert password.")
    system("rm -f ~/edm/weight-pred/model/run/*")
    system( "scp -r daon@access.cims.nyu.edu:~/weight-pred/model/runs/ ." )
}

load( "runs/parameters.Rdata" )

## Get the true values we were trying to predict
df <- read.csv("originals/three_species.csv",
                   header = TRUE,
                   sep = "," )
df <- lag_every_variable(df, n_lags)
truth <- df$x_p1[pred[1]:pred[2]]

weighted_pred <- read.table("runs/weighted_predictions_x.csv",
                            na = "NA",
                            sep =",")
                            ## row.names = FALSE,
                            ## col.names = FALSE)

mve_pred <- read.table("runs/mve_predictions_x.csv",
                       na = "NA",
                       sep ="," )

## Calculate mean prediction skill per library size
weighted_rho <- mean_cor( weighted_rho, truth, nrow = n_samp) )
mve_rho      <- mean_cor(      mve_rho, truth, nrow = n_samp) )

print( weighted_rho )
print( mve_rho )
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
