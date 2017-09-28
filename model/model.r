#!/usr/bin/Rscript
library(rEDM)
library(zoo)
source( "../helpers/helper.r" )

## Get a weighted prediction predictors and their uncertainty
weighted_prediction <- function(var_table,
                                pred_table)
{
    ## Sanity check
    stopifnot(
        nrow(var_table) == nrow(pred_table) &&
        ncol(var_table) == ncol(pred_table)
    )

    ## Meomry allocation
    ret <- numeric(ncol(pred_table))

    ## Avoid loops at all costs!!! Or be lazy!!
    for( i in 1:ncol(pred_table) )
    {            
        ## Extract uncertainty of all predictors and find its
        v   <- var_table[ , i]
        ind <- order( v )[1:ceiling(sqrt(nrow(var_table)))]
        
        ## Find who is in, get corresponding predictions and get
        ## corresponding weights, both exponential and precision.
        p <- pred_table[ind, i]
        w <- exp(-v[ind])
            
        ## Weight the predictors according to the weights and store
        ## the weighted prediction in the preallocated matrices.
        ret[i] <- sum(p * w) / sum(w)
    }
    return( ret )
}


## If you want to download the file for some reason... don't!
## filename = paste0("http://science.sciencemag.org/highwire/filestream/683325/",
##                   "field_highwire_adjunct_files/1/aag0863_SupportingFile_Suppl._Excel_seq1_v2.xlsx")

## Load data and immediately normalize to zero mean and unit std dev
df <-normalize_df(read.csv("originals/data.csv",
                           header = TRUE,
                           sep = "," )
                  )

## Parameters
n_vars <- ncol(df) ## == 3 cuz x, y, z
E <- 3
n_lags <- 3 ## Lags: 0, -1, -2
n_samp <- 3 ##100 ## Number of random libraries
pred <- c(2501,3000) ## Prediciton always constant
pred_size <- pred[2] - pred[1] + 1
lib_sizes <- (1:3)*10

## Library changes in size and is also random
results <- data.frame(lib_sizes = lib_sizes)
results$avg <- numeric(nrow(results))
results$top <- numeric(nrow(results))
results$med <- numeric(nrow(results))
results$bot <- numeric(nrow(results))

## Craeate lags for every variable and make y_p1
df <- lag_every_variable(df, n_lags)
y_p1 <- c( as.vector(lag(zoo(df$y),1)) , NA )
df <- data.frame(y_p1 = y_p1, df )


## We know the allowed combinations must have one unlagged
## variable. By the way function combn is structured, these will be
## the first n_comb (defined below).
combinations <- combn( n_vars*n_lags, E ) ## Get ALL cobminations.
## Only this many are OK.
n_comb <- choose( n_vars*n_lags, E )- choose(n_vars*(n_lags-1), E ) 
## Throw away the rest.
combinations <- combinations[ ,1:n_comb ]

## Preallocate memory for predictions and uncertainties
pred_table  <- matrix( NA,  nrow = n_comb, ncol = pred_size )
var_table   <- matrix( Inf, nrow = n_comb, ncol = pred_size )
predictions <- matrix( NA,  nrow = n_samp, ncol = pred_size )
rhos        <- rep   ( NA,         n_samp                  )



for (lib_ind in 1:nrow(results) )
{
    ## For every random library starting point
    for( smp in 1:n_samp )
    {
        ## Choose random lib of length exactly 100. Pred is always the same
        ## and they never overlap.
        lib  <- sample(501:2001, 1)
        lib  <- c(lib, lib + results$lib_size[lib_ind]- 1 )
        
        for( i in 1:n_comb )
        {
            cols <- combinations[,i]
            
            ## Just to make sure, all cols have one unlagged coordinate.
            stopifnot( min(cols) <= n_vars )
            cols <- cols + 1 ## cuz first col is y_p1
            
            output <- block_lnlp(df,
                                 lib = lib,   ## lib is chosen randomly above
                                 pred = pred, ## pred is ALWAYS the same
                                 method = "simplex",
                                 tp = 0, ## Zero step cuz it is built into y_p1
                                 columns = cols, 
                                 target_column = 1, ## == y_p1
                                 first_column_time = FALSE,
                                 short_output = TRUE,
                                 stats_only = FALSE )

            time <- output[[1]]$model_output$time
            pr   <- output[[1]]$model_output$pred
            vars <- output[[1]]$model_output$pred_var
            rho <- output[[1]]$stats$rho
                    
            calc_rho <- cor(pr, y_p1[time], use = "complete.obs" ) 
            
            ## Not sure why these are not exactly equal...
            ## if(  abs(calc_rho-rho) > 1e-14 )
            ##     stop( calc_rho - rho )

            ## Save ALL data.
            pred_table[i, ] <- pr 
            var_table [i, ] <- vars
            
        }
        
        predictions[smp, ] <- weighted_prediction(var_table, pred_table)
        rhos[smp] <- cor(predictions[smp, ], y_p1[time], use = "complete.obs" )
    }

    quarts <- quantile(rhos, probs = c(0.25,0.5,0.75))

    results$avg[lib_ind] <- mean(rhos)
    results$bot[lib_ind] <- quarts[1]
    results$med[lib_ind] <- quarts[2]
    results$top[lib_ind] <- quarts[3]
}

write.csv(results,
          file = "data/rhos.csv",
          row.names = FALSE,
          col.names = TRUE)


## noo <- read.csv("data/rhos.csv", header = TRUE, sep = "," )
## print(noo)
