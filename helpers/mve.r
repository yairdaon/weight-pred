#!/usr/bin/Rscript
library(rEDM)
library(zoo)
source( "helpers/helper.r" )

## Get a uncertainty-weighted prediction Variable num_neighbors is
## redundant, but we keep for consistency of method signature. Also
## k is currently redundant.
uwe <- function(df,
                lib = c(1, floor(NROW(block)/2)), 
                pred = c(floor(NROW(block)/2), NROW(block)),
                norm_type = c("L2 norm", "L1 norm", "P norm"),
                P = 0.5, 
                E = 3,
                tau = 1,
                tp = 1,
                max_lag = 3,
                num_neighbors = E+1, 
                k = "sqrt",
                na.rm = FALSE,
                target_column = 1, 
                stats_only = TRUE,
                first_column_time = FALSE,
                exclusion_radius = NULL,
                silent = FALSE, 
                short_output = FALSE
                )
                
{
    ## combinations <- best_combinations_cv(df,
    ##                                      lib,
    ##                                      target_column,
    ##                                      max_lag,
    ##                                      E)
    combinations <- name_combinations(df,max_lag,E)
        
    df <- lag_every_variable(df, max_lag)
    df <- respect_lib( df, lib, max_lag )
    
    ## Useful variables to have
    n_comb  <- ncol(combinations)
    pred_size <- pred[2] - pred[1]
    
    ## Preallocate memory for everything.
    pred_table <- matrix( NA,  nrow = n_comb, ncol = pred_size )
    var_table  <- matrix( Inf, nrow = n_comb, ncol = pred_size )
        
    for( comb in 1:n_comb )
    {
        output <- block_lnlp(df,
                             lib  = lib,  
                             pred = pred, 
                             method = "simplex",
                             tp = 1, 
                             columns = combinations[ , comb ], 
                             num_neighbors = E+1, ## Do NOT use 1!!
                             target_column = target_column, 
                             first_column_time = FALSE,
                             stats_only = FALSE )
        pred_table[ comb, ] <- output$model_output[[1]]$pred
        var_table [ comb, ] <- output$model_output[[1]]$pred_var
        
    } ## Closes for( comb in 1:n_comb )
    
    ## Get rid of those zero variance predictions.
    bad_var <- var_table == 0
    var_table[ bad_var ] <- Inf
    
    ## Apply order to every column. Use increasing order, so lower
    ## uncertainties have lower indices.
    ord <- colOrder( var_table )

    ## Now var_table and pred_table both have columns ordered
    ## according to var_table.
    var_table  <- matrix(var_table [ord], ncol = ncol(var_table))
    pred_table <- matrix(pred_table[ord], ncol = ncol(var_table))

    ## Take best sqrt(k), according to the MVE heuristic
    best <- ceiling(sqrt(nrow(var_table)))

    ## Take only best predictors
    var_table  <- var_table [1:best, ]
    pred_table <- pred_table[1:best, ]
    
    ## Exponential weights
    weight_table <- exp(-var_table )
    
    ## Weighted sum
    predictions <- colSums(weight_table*pred_table) / colSums(weight_table)
    time        <- output$model_output[[1]]$time
    obs         <- output$model_output[[1]]$obs

    params <- data.frame(E = E, tau = tau, tp = tp, nn = num_neighbors) ##, k = k_list)
    stats <- compute_stats(obs, predictions)

    if( stats_only )
        return( cbind(params, stats) )
    else 
        return(list(params = params, 
                    model_output = data.frame(time = time, 
                                              obs = obs, 
                                              pred = predictions), 
                    pred_stats = stats )
               )
    
}

mve <- function(df,
                lib = c(1, floor(NROW(block)/2)), 
                pred = c(floor(NROW(block)/2), NROW(block)),
                norm_type = c("L2 norm", "L1 norm", "P norm"),
                P = 0.5, 
                E = 3,
                tau = 1,
                tp = 1,
                max_lag = 3,
                num_neighbors = E+1,
                k = "sqrt",
                na.rm = FALSE,
                target_column = 1,
                stats_only = TRUE,
                first_column_time = FALSE,
                exclusion_radius = NULL,
                silent = FALSE, 
                short_output = FALSE
                )
    
{
    ## Get the best models to be used by MVE.
    best <- best_combinations_cv(df,
                                 lib,
                                 target_column,
                                 max_lag,
                                 E)
        
    ## Check: respects lib? Order of vars.
    df <- lag_every_variable(df, max_lag)
    df <- respect_lib( df, lib, max_lag )
    
    predictions <- 0
    
    ## Find prediction for every one of the best models
    for( comb in 1:ncol(best) )
    {
        cols <- best[ , comb ]
        output <- block_lnlp(df,
                             lib  = lib,
                             pred = pred,
                             method = "simplex",
                             tp = tp,
                             columns = cols,
                             num_neighbors = num_neighbors,
                             target_column = target_column, 
                             first_column_time = FALSE,##first_column_time,
                             stats_only = FALSE,
                             silent = FALSE,
                             exclusion_radius = exclusion_radius )
        predictions <- predictions + output$model_output[[1]]$pred
    }
    
    ## The predictions MVE makes
    predictions <- predictions / ncol(best)
    obs  <- output$model_output[[1]]$obs
    time <- output$model_output[[1]]$time ## Not sure what this does...
    
    params <- data.frame(E = E, tau = tau, tp = tp, nn = num_neighbors) ##, k = k_list)
    stats <- compute_stats(obs, predictions)

    if( stats_only )
        return( cbind(params, stats) )
    else 
        return(list(params = params, 
                    model_output = data.frame(time = time, 
                                              obs  = obs, 
                                              pred = predictions), 
                    pred_stats = stats )
               )
    
}

## Finds the best models according to the MVE heuristic
best_combinations_cv <- function(df,
                                 lib,
                                 target_column,
                                 max_lag,
                                 E)
{
    combinations <- name_combinations(df,max_lag,E)
    n_comb <- ncol(combinations)
   
    df <- lag_every_variable(df, max_lag)
    df <- respect_lib( df, lib, max_lag )
    
    ## Preallocate memory.
    rhos <- numeric( n_comb )
    for( comb in 1:n_comb )
    {
        cols <- combinations[ , comb ]
        output <- block_lnlp(df,
                             lib  = lib,
                             pred = lib, ## CV!!!
                             method = "simplex",
                             tp = 1,
                             columns = cols,
                             target_column = target_column,
                             first_column_time = FALSE,
                             stats_only = TRUE,
                             silent = TRUE )
        rhos[comb] <- output$rho
    }
    
    ord <- order( rhos, decreasing = TRUE )
    ind <- ord[ 1 : ceiling(sqrt(length(rhos))) ]
    ret <- combinations[ , ind ]
    
    return( ret )
}

