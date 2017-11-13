#!/usr/bin/Rscript
library(rEDM)
library(zoo)
source( "../helpers/helper.r" )


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
##                                 short_output = FALSE,
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
                method = "uniform")       
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
                             short_output = FALSE,
                             stats_only = FALSE )
        pred_table[comb, ] <- output[[1]]$model_output$pred    [pred[1]:pred[2]] ## Prediction 
        var_table [comb, ] <- output[[1]]$model_output$pred_var[pred[1]:pred[2]] ## Uncertainty

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
    
    prediction <- colSums(weights * pred_table) / colSums(weights) ## Weighted average

    ## Debugging...
    if( any( is.na(prediction) ) ) {
        print( paste0( "NAs:", sum(is.na(prediction)) ) )
        print( "pred table:" )
        print( pred_table )
        print( "var table:" )
        print( var_table )
        print( "pred:" )
        print( prediction )
        print( lib )
        print("")
    }
    return( prediction )
           
} ## Closes we <- function(...)

