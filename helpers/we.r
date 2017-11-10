#!/usr/bin/Rscript
library(rEDM)
library(zoo)
source( "helper.r" )

weighted_colMeans <- function(X,
                              V ) ## variances / uncertainties
                              
{
    W <- 1 / V ## Unnormalized weights 
    return( colSums(W * X) / colSums(W) ) ## Weighted average
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
    P <- matrix( NA,  nrow = n_comb, ncol = pred_size )
    V <- matrix( NA,  nrow = n_comb, ncol = pred_size )
            
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
        P[comb, ] <- output[[1]]$model_output$pred    [pred[1]:pred[2]] ## Prediction 
        V[comb, ] <- output[[1]]$model_output$pred_var[pred[1]:pred[2]] ## Uncertainty
                
    } ## Closes for( comb in 1:n_comb )

    prediction <- weighted_colMeans(P,V)
    return( prediction )
           
} ## Closes we <- function(...)

