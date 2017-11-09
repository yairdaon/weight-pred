#!/usr/bin/Rscript
library(rEDM)
library(zoo)
source( "../helpers/helper.r" )

weighted_colMeans <- function(X,
                              W, ## Weights, NOT variances / uncertainties
                              n = ceiling(sqrt(nrow(W))) )
{
    ret <- numeric(ncol(X))
        
    ## Avoid loops at all costs!!! Or be lazy!!
    for( i in 1:ncol(X) )
    {            
        w   <- W[ , i]
        ind <- order( w, decreasing = TRUE )[1:n]
        w   <- w[ind]
        
        p <- X[ind, i]
                    
        ret[i] <- sum(p * w) / sum(w)
    }

    return(ret)
}

we <- function(df, ## with lagged and scaled variables.
               variable, ## For which we wish to make predictions.
               lib,  ## Library set.
               pred, ## Prediciton set.
               combinations ## Produced by make_combintaions.
               )       
{
    ## Useful variables to have
    target  <- paste0( variable, "_p1" )
    n_comb  <- ncol(combinations)
    pred_size <- pred[2] - pred[1] + 1
    
    ## Preallocate memory for everything.
    P    <- matrix( NA,  nrow = n_comb, ncol = pred_size )
    W    <- matrix( NA,  nrow = n_comb, ncol = pred_size )
            
    for( comb in 1:n_comb )
    {
        ## Variables used for prediction are chosen
        ## according to the current combinations. Since
        ## first cols are x_p1, y_p1, z_p1, we skip them
        ## and add n_var.
        cols   <- combinations[ , comb ] ## see ***
        
        output <- block_lnlp(df,
                             lib  = lib,  
                             pred = pred, 
                             method = "simplex",
                             tp = 0, ## Lags built into df
                             columns = cols, 
                             target_column = target, 
                             first_column_time = FALSE,
                             short_output = FALSE,
                             stats_only = FALSE )
        P[comb, ] <-      output[[1]]$model_output$pred    [pred[1]:pred[2]] 
        W[comb, ] <- exp(-output[[1]]$model_output$pred_var[pred[1]:pred[2]])
                
    } ## Closes for( comb in 1:n_comb )

    prediction <- weighted_colMeans(P,W)
    return( prediction )
           
} ## Closes we <- function(...)

