#!/usr/bin/Rscript
library(rEDM)
library(zoo)

## If you want to download the file for some reason... don't!
## filename = paste0("http://science.sciencemag.org/highwire/filestream/683325/",
##                   "field_highwire_adjunct_files/1/aag0863_SupportingFile_Suppl._Excel_seq1_v2.xlsx")


df <- read.csv("originals/data.csv",
               header = TRUE,
               sep = "," )

## Lags: 0, -1, -2
n_lags <- 3 

## Number of random libraries
n_samp <- 2

## Craeate lags for every variable
for( var in names(df) )
{
    ts <- zoo( df[ , var ] )
    for( k in 1:(n_lags-1) )
        df[ paste0(var, "_", k) ] <- c(rep( NA, k ),
                                       lag( ts, -k )
                                       )
}

## Normalize to zero mean and unit std dev
vars <- names(df) 
df <- scale(df) ## This makes the df a matrix...
df <- data.frame( df ) ## ...so make it a df...
names(df) <- vars ##...and restore the names.


## We know the allowed combinations must have one unlagged
## variable. By the way function combn is structured, these will be
## the first n_comb (defined below).
combinations <- combn( length(vars), 3 ) ## Get ALL cobminations.
n_comb <- choose(3*n_lags,3) - choose(3*(n_lags-1),3) ## Only this
                                                      ## many are OK.
combinations <- combinations[ ,1:n_comb ] ## Throw away the rest.

## Preallocate memory for predictions and uncertainties
pred_table <- matrix( NA,  nrow = n_comb*n_samp, ncol = 500 )
var_table  <- matrix( Inf, nrow = n_comb*n_samp, ncol = 500 )

## For every random library starting point
for( smp in 0:(n_samp-1) )
{
    ## Choose random lib of length exactly 100. Pred is always the same
    ## and they never overlap.
    lib  <- sample(501:2001, 1)
    lib  <- c(lib, lib + 99 )
    
    
    
    for( i in 1:n_comb )
    {
        cols <- combinations[,i]
        
        ## Just to make sure, all cols have one unlagged coordinate.
        stopifnot( min(cols) <= 3 )
        
        output <- block_lnlp(df,
                             lib = lib,
                             pred = c(2501,3000), ## pred is ALWAYS the same
                             method = "simplex",
                             tp = 1, ## Exactly one step ahead
                             columns = cols, 
                             target_column = 2, ## the y coordinate
                             first_column_time = FALSE,
                             short_output = TRUE,
                             stats_only = FALSE )
        
        ## Get model output: prediction times, predictions and variances:
        time <- output[[1]]$model_output$time - pred[1] + 1
        pr   <- output[[1]]$model_output$pred
        vars <- output[[1]]$model_output$pred_var
        
        ## Sanity checks #####################
        len_time <- length(time)
        len_vars <- length(vars)
        len_pred <- length(pr)
        
        stopifnot( len_time == len_vars )
        stopifnot( len_time == len_pred )

        no_vars <- is.na(vars) | is.nan(vars)
        no_pr   <- is.na(pr) 
        if( any(no_vars!=no_pr) )
            stop( "No prediction indices do not align with no variance indices!" )
        ind <- no_pr
        
        pr[ ind ] <- NA
        vars[ ind ] <- Inf

        ## Save ALL data.
        pred_table[ i + n_comb*smp, time ] <- pr
        var_table [ i + n_comb*smp, time ] <- vars    
    }
}
