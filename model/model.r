#!/usr/bin/Rscript
library(rEDM)
library(zoo)
source( "../helpers/helper.r" )


## In this script we show prediction skill of the uncertainty -
## weighting scheme. We randomly sample libraries from synthetic data
## and present their skill at out-of-sample prediction.

## If you want to download the file for some reason... don't!
## filename = paste0("http://science.sciencemag.org/highwire/filestream/683325/",
##                   "field_highwire_adjunct_files/1/aag0863_SupportingFile_Suppl._Excel_seq1_v2.xlsx")

save_weighted_predictions <- function(raw_df,
                                      variables = names(raw_df), 
                                      E = 3, ## Embedding dimension of the system.
                                      n_lags = 3, ## 0, -1,..., -(n_lags-1)
                                      n_samp = 5, ## Number of random libraries, should be in the hundreds
                                      lib_range = c(501:2001),  ## Library set.
                                      pred = c(2501,3000), ## Prediciton set.
                                      lib_sizes = (1:3)*10     ## Library changes in size and is also random
                                      )
{   
    ## Immediately rescale the data frame and keep track of the
    ## normalizing factors.
    mus  <- get_mus(raw_df)
    sigs <- get_sigs(raw_df)
    df   <- data.frame(scale(raw_df))
      
    ## number of variables in the dynamical system
    n_vars <- ncol(df)

    ## Create lags for every variable and make y_p1
    df <- lag_every_variable(df, n_lags)

    ## Create the list of combinations
    combinations <- make_combinations(n_vars, n_lags, E)
    n_comb       <- ncol(combinations)
    
    ## Preallocate memory for predictions and uncertainties
    pred_size   <- pred[2] - pred[1] + 1
    pred_table  <- matrix( NA,  nrow = n_comb, ncol = pred_size )
    var_table   <- matrix( Inf, nrow = n_comb, ncol = pred_size )
    predictions <- matrix( NA,  nrow = n_samp, ncol = pred_size )

    ## Do the analysis for every variable. 
    for( curr_var in variables )
    {
        for (lib_size in lib_sizes)
        {
            ## For every random library starting point
            for( smp in 1:n_samp )
            {
                ## Choose random lib
                lib <- random_lib(lib_range, lib_size)
                
                for( comb in 1:n_comb )
                {

                    ## Just to make sure, all variable sets have one
                    ## unlagged coordinate.
                    ## stopifnot( min(combinations[,comb]) <= n_vars )
                    
                    ## Variables used for prediction are chosen
                    ## according to the current combinations. Since
                    ## first cols are x_p1, y_p1, z_p1, we skip them
                    ## and add n_var.
                    output <- block_lnlp(df,
                                         lib = lib,   ## lib is chosen randomly above
                                         pred = pred, ## pred is ALWAYS the same
                                         method = "simplex",
                                         tp = 0, ## Lags built into df
                                         columns = combinations[, comb] + n_vars, 
                                         target_column = paste0( curr_var, "_p1" ),
                                         first_column_time = FALSE,
                                         short_output = FALSE,
                                         stats_only = FALSE )

                    ## time <- output[[1]]$model_output$time
                    ## pr   <- output[[1]]$model_output$pred[pred[1]:pred[2]]
                    ## vars <- output[[1]]$model_output$pred_var[pred[1]:pred[2]]

                    ## rho <- output[[1]]$stats$rho
                    ## calc_rho <- cor(pr, df[time,var_ind], use = "complete.obs" ) 
                    ## Not sure why these are equal only up to 1e-15 error...
                    ## if(  abs(calc_rho-rho) > 1e-14 )
                    ##     stop( calc_rho - rho )

                    ## These tables hold data for the current (random) library
                    pred_table[comb, ] <- output[[1]]$model_output$pred    [pred[1]:pred[2]] 
                    var_table [comb, ] <- output[[1]]$model_output$pred_var[pred[1]:pred[2]]
                }

                ## Weighted predictions for the current sampled library 
                predictions[smp, ] <- weighted_prediction(var_table, pred_table)
                ## prediction <- descale(prediction, mus$curr_var, sigs$curr_var)
            }

            write.table(predictions,
                        file = paste0("data/all_predictions_", curr_var, "_lib_size_", lib_size, ".csv"),
                        quote = FALSE,
                        na = "NA",
                        row.names = FALSE,
                        col.names = FALSE)
        }
    }
}


## Load data 
df <- read.csv("originals/data.csv",
                     header = TRUE,
                     sep = "," )

variables <- c( "x", "y", "z" )

save_weighted_predictions(df,
                          variables = variables,
                          E = 3, ## Embedding dimension of the system.
                          n_lags = 3, ## 0, -1,..., -(n_lags-1)
                          n_samp = 150, ## Number of random libraries, should be in the hundreds
                          lib = c(501,2001),  ## Library set.
                          pred = c(2501,3000), ## Prediciton set.
                          lib_sizes = (1:30)*5 ## Library sizes
                          )
