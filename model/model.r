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

predictions <- function(filename = stop("File name must be provided!"),
                        variables = names(raw_df), 
                        E = 3, ## Embedding dimension of the system.
                        n_lags = 3, ## 0, -1,..., -(n_lags-1)
                        n_samp = 200, ## Number of random libraries, should be in the hundreds
                        lib = c(501:2001),  ## Library set.
                        pred = c(2501,3000), ## Prediciton set.
                        lib_sizes = (1:15)*10     ## Library changes in size and is also random
                        )
{
    ## Load data 
    raw_df <- read.csv(filename,
                   header = TRUE,
                   sep = "," )

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
    
    ## Preallocate memory for everything.
    pred_size   <- pred[2] - pred[1] + 1
    pred_table  <- matrix( NA,  nrow = n_comb,                   ncol = pred_size )
    var_table   <- matrix( Inf, nrow = n_comb,                   ncol = pred_size )
    predictions <- matrix( NA,  nrow = n_samp*length(lib_sizes), ncol = pred_size )
    mve         <- matrix( NA,  nrow = n_samp*length(lib_sizes), ncol = pred_size )
    rand_libs   <- matrix( NA,  nrow = n_samp*length(lib_sizes), ncol = 2         )
    rhos        <- numeric( n_comb )
    
    ## Do the analysis for every variable. 
    for( curr_var in variables )
    {
        lib_ind = 0
        for (lib_size in lib_sizes)
        {

            ## For every random library starting point
            for( smp in 1:n_samp )
            {
                ## Choose random lib
                rand_lib <- random_lib(lib, lib_size)
                lib_ind <- lib_ind + 1
                rand_libs[ lib_ind, ] <- rand_lib
                    
                for( comb in 1:n_comb )
                {

                    ## Just to make sure, all variable sets have one
                    ## unlagged coordinate.
                    ## stopifnot( min(combinations[,comb]) <= n_vars )
                    
                    ## *** Variables used for prediction are chosen
                    ## according to the current combinations. Since
                    ## first cols are x_p1, y_p1, z_p1, we skip them
                    ## and add n_var.
                    output <- block_lnlp(df,
                                         lib = rand_lib, ## chosen randomly above
                                         pred = pred, ## always the same
                                         method = "simplex",
                                         tp = 0, ## Lags built into df
                                         columns = combinations[, comb] + n_vars, ## see ***
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


                    ## Now we do CV, for the MVE ranking. All we need
                    ## is the rho value for ranking.
                    cv_output <- block_lnlp(df,
                                         lib = rand_lib, ## chosen randomly above
                                         pred = rand_lib, ## CV!!!
                                         method = "simplex",
                                         tp = 0, ## Lags built into df
                                         columns = combinations[, comb] + n_vars, ## see ***
                                         target_column = paste0( curr_var, "_p1" ),
                                         first_column_time = FALSE,
                                         stats_only = TRUE )
                    rhos[comb] <- cv_output$rho
                    
                } ## Closes  for( comb in 1:n_comb )

                ## Weighted predictions for the current sampled library 
                prediction <- weighted_prediction(var_table, pred_table)
                prediction <- descale(prediction, mus$curr_var, sigs$curr_var)
                predictions[lib_ind, ] <- prediction

                
                ## MVE Predictions:
                ord <- order(rhos)[1:ceiling(sqrt(n_comb))]
                mve[ lib_ind, ] <- colMeans( pred_table[ord, ], na.rm = TRUE ) 
                
            } ## Closes for( smp in 1:n_samp )
            
        } ## Closes for (lib_size in lib_sizes)


        write.table(predictions,
                    file = paste0("runs/weighted_predictions_", curr_var, ".csv"),
                    quote = FALSE,
                    na = "NA",
                    row.names = FALSE,
                    col.names = FALSE)

        write.table(mve,
                    file = paste0("runs/mve_predictions_", curr_var, ".csv"),
                    quote = FALSE,
                    na = "NA",
                    row.names = FALSE,
                    col.names = FALSE)
        
        write.table(rand_libs,
                    file = paste0("runs/libraries_", curr_var, ".txt" ),
                    quote = FALSE,
                    na = "NA",
                    row.names = FALSE,
                    col.names = FALSE,
                    sep = "\t")
                    
        
    } ## Closes for( curr_var in variables )
    
    ## Now that we are done, we save the parameters so we can know
    ## EXACTLY what parameters we were using
    save(filename,
         variables,
         E,
         n_lags,
         n_samp,
         lib,
         pred,
         lib_sizes,         
         file = "runs/parameters.Rdata" )

} ## Closes predictions <- function(...)



predictions(file = "originals/three_species.csv",
            variables =c( "x", "y" ),
            E = 3, ## Embedding dimension of the system.
            n_lags = 3, ## 0, -1,..., -(n_lags-1)
            n_samp = 5, ## Number of random libraries, should be in the hundreds
            lib = c(501,2001),  ## Library set.
            pred = c(2501,3000), ## Prediciton set.
            lib_sizes = (1:2)*25 ## Library sizes
            )


