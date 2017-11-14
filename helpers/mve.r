#!/usr/bin/Rscript
library(rEDM)
library(zoo)
source( "helpers/helper.r" )


## find_best_fraction <- function(df,
##                                variable,
##                                combinations,
##                                lib)
## {
##     n_comb     <- ncol(combinations)
##     pred_table <- matrix( NA,  nrow = n_comb, ncol = lib[2] - lib[1] + 1 )
##     var_table  <- matrix( NA,  nrow = n_comb, ncol = lib[2] - lib[1] + 1 )
##     target     <- paste0( variable, "_p1" )
    
##     for( comb in 1:ncomb )
##     {
##         cols   <- combinations[ , comb ] ## see ***
##         cv_output <- block_lnlp(df,
##                                 lib  = lib, 
##                                 pred = lib, ## CV!!!
##                                 method = "simplex",
##                                 tp = 0, ## Lags built into df
##                                 columns = cols,
##                                 target_column = target,
##                                 first_column_time = FALSE,
##                                 stats_only = FALSE,
##                                 silent = TRUE )
##         pred_table[comb, ] <- output[[1]]$model_output$pred    [lib[1]:lib[2]] ## Prediction 
##         var_table [comb, ] <- output[[1]]$model_output$pred_var[lib[1]:lib[2]] ## Uncertainty
##     }

##     pred <- 
## }




mve <- function(df, ## with lagged and scaled variables.
                variable, ## For which we wish to make predictions.
                lib,  ## Library set.
                pred, ## Prediciton set.
                combinations, ## Produced by make_combintaions.
                method)       
{
    ## Useful variables to have
    target  <- paste0( variable, "_p1" )
    n_comb  <- ncol(combinations)
    pred_size <- pred[2] - pred[1] + 1
    
    ## Preallocate memory for everything.
    pred_table <- matrix( NA,  nrow = n_comb, ncol = pred_size )
    var_table  <- matrix( Inf, nrow = n_comb, ncol = pred_size )
    rhos       <- numeric( n_comb )
    
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
                             stats_only = FALSE )
        ## times <- output$model_output[[1]]$time
        ## troof <- output$model_output[[1]]$obs
        
        pred_table[ comb, ] <- output$model_output[[1]]$pred
        var_table [ comb, ] <- output$model_output[[1]]$pred_var

        ## bad <- which( var_table[ comb, ]  == 0  )
        ## if( length( bad ) > 0 )
        ## {
        ##     print( lib )
        ##     print( bad )
        ## }

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

    ord  <- order( rhos, decreasing = TRUE )
    best <- ord[ 1 : ceiling(sqrt(length(ord))) ]
    pred_table <- pred_table[ best, ]
    var_table  <- var_table [ best, ]
    
    if( method == "uniform" )
        weights <- matrix( 1, ncol = ncol(var_table), nrow = nrow(var_table) )
    else if ( method == "minvar" )
        weights <- 1 / var_table 
    else if ( method == "exp" )
        weights <- exp( -var_table )
    else
        stop( paste0("Unknown weighting scheme '", method, "'" ) )

    ## Fix bad stuff. This is ad-hoc and I don't like it.
    bad <- which( (var_table == 0) || is.na(pred_table) || is.nan(pred_table) )
    weights   [ bad ] <- 0
    pred_table[ bad ] <- 0
    
    prediction <- colSums(weights * pred_table) / colSums(weights) ## Weighted average
    return( prediction )
           
} ## Closes we <- function(...)

