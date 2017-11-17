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

## Get a uncertainty-weighted prediction
uwe <- function(df, ## lagged and scaled.
                variable,
                lib, ## Library set.
                pred, ## Prediciton set.
                combinations)

{
    ## Useful variables to have
    target  <- paste0( variable, "_p1" )
    n_comb  <- ncol(combinations)
    pred_size <- pred[2] - pred[1] + 1
    
    ## Preallocate memory for everything.
    pred_table <- matrix( NA,  nrow = n_comb, ncol = pred_size )
    var_table  <- matrix( Inf, nrow = n_comb, ncol = pred_size )
        
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
        pred_table[ comb, ] <- output$model_output[[1]]$pred
        var_table [ comb, ] <- output$model_output[[1]]$pred_var
        
    } ## Closes for( comb in 1:n_comb )

    truth <- output$model_output[[1]]$obs
    
    ## Get rid of those zero variance predictions.
    var_table[ var_table == 0 ] <- NA
    
    ## Get increasing order, so lower uncertainties
    ## have lower indices.
    ord <- colOrder( var_table )

    ## Take best sqrt(k), according to the MVE heuristics
    best <- ceiling(sqrt(nrow(var_table)))
    ind  <- ord[1:best, ]

    w <- exp(-matrix(var_table [ind], nrow = best, ncol = ncol(var_table) ) )
    p <-      matrix(pred_table[ind], nrow = best, ncol = ncol(var_table) )
    ret <- colSums( w*p ) / colSums(w)

    attr(ret, "rho") <- cor( ret, truth, use = "complete.obs" )
    ## Set attribute to ret of rho
    return( ret )
}

mve <- function(df, ## with lagged and scaled variables.
                variable, ## For which we wish to make predictions.
                lib,  ## Library set.
                pred, ## Prediciton set.
                combinations) ## Produced by make_combintaions.       
{
    ## Useful variables to have
    target  <- paste0( variable, "_p1" )
    n_comb  <- ncol(combinations)
    pred_size <- pred[2] - pred[1] + 1
    
    ## Preallocate memory for everything.
    pred_table <- matrix( NA,  nrow = n_comb, ncol = pred_size )
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
        pred_table[ comb, ] <- output$model_output[[1]]$pred
        
        rhos[comb] <- block_lnlp(df,
                                lib  = lib, 
                                pred = lib, ## CV!!!
                                method = "simplex",
                                tp = 0, ## Lags built into df
                                columns = cols,
                                target_column = target,
                                first_column_time = FALSE,
                                stats_only = TRUE,
                                silent = TRUE )$rho
        
    } ## Closes for( comb in 1:n_comb )

    truth <- output$model_output[[1]]$obs
        
    ord <- order( rhos, decreasing = TRUE )
    ind <- ord[ 1 : ceiling(sqrt(length(rhos))) ]
    ret <- colMeans( pred_table[ ind, ] )

    attr(ret, "rho") <- cor( ret, truth, use = "complete.obs" )
    ## Set attribute to ret of rho
    return( ret )
}

