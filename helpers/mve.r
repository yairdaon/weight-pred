#!/usr/bin/Rscript
library(rEDM)
library(zoo)

mve <- function(df, ## with lagged and scaled variables.
                variable, ## For which we wish to make predictions.
                lib,  ## Library set.
                pred, ## Prediciton set.
                combinations ## Produced by make_combintaions.
                )       
{
    ## Useful variables to have
    target  <- paste0( variable, "_p1" )
    n_comb  <- ncol(combinations)
    
    ## Preallocate memory for everything.
    pred_table  <- matrix( NA,  nrow = n_comb, ncol = pred[2] - pred[1] + 1 )
    rhos        <- numeric( n_comb )
        
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
        pred_table[comb, ] <- output[[1]]$model_output$pred[pred[1]:pred[2]] 
        
        ## Now we do CV, for the MVE ranking. All we need
        ## is the rho value for ranking.
        cv_output <- block_lnlp(df,
                                lib  = lib, 
                                pred = lib, ## CV!!!
                                method = "simplex",
                                tp = 0, ## Lags built into df
                                columns = cols,
                                target_column = target,
                                first_column_time = FALSE,
                                stats_only = TRUE,
                                silent = TRUE )
        rhos[comb] <- cv_output$rho
        
    } ## Closes for( comb in 1:n_comb )
    
    ord <- order( rhos, decreasing = TRUE )
    top_performers <- pred_table[ord[ 1 : ceiling(sqrt(length(rhos))) ] , ]
    prediction <- colMeans(top_performers, na.rm = TRUE ) 
    return( prediction )
           
} ## Closes mve <- function(...)

