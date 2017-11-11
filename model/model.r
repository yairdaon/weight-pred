#!/usr/bin/Rscript
library(rEDM)
library(zoo)
source( "../helpers/helper.r" )

save_predictions <- function(filename = stop("File name must be provided!"),
                             variables = names(raw_df), 
                             E = 3, ## Embedding dimension of the system.
                             n_lags = E, ## 0,-1, ..., -n_lags
                             n_samp = 200, ## Number of random libraries, should be in the hundreds
                             lib = c(501:2001),  ## Library set.
                             pred = c(2501,3000), ## Prediciton set.
                             lib_sizes = (1:15)*10,     ## Library changes in size and is also random
                             method = "we"
                             )
{
    if( method == "we" )
    {
        source( "../helpers/we.r" )
        pred_func <- we
    }
    else if( method == "mve" )
    {
        source( "../helpers/mve.r" )
        pred_func <- mve
    }
    else
        stop( "Choose method to be 'me' or 'mve'" )
    
    ## Load data 
    raw_df <- read.csv(filename,
                       header = TRUE,
                       sep = "," )
    
    ## Rescale the data frame and keep track of the normalizing
    ## factors.
    mus  <- get_mus(raw_df)
    sigs <- get_sigs(raw_df)
    df   <- data.frame(scale(raw_df))    
    
    ## Create lags for every variable and make y_p1
    df <- lag_every_variable(df, n_lags)
        
    ## Create the list of combinations
    combinations <- make_combinations(ncol(raw_df), ## total number of variables
                                      n_lags,
                                      E)
            
    ## Preallocate memory for everything.
    predictions <- matrix(NA,
                          nrow = n_samp*length(lib_sizes),
                          ncol = pred[2] - pred[1] + 1 )
    
    ## Do the analysis for every variable. 
    for( curr_var in variables )
    {
        lib_ind <- 0
        for (lib_size in lib_sizes)
        {

            ## For every random library starting point
            for( smp in 1:n_samp )
            {
                ## Choose random lib
                rand_lib <- random_lib(lib, lib_size)
                lib_ind  <- lib_ind + 1    

                ## Find the MVE prediction
                prediction <- pred_func(df, ## lagged and scaled.
                                        curr_var,
                                        rand_lib, ## Library set.
                                        pred, ## Prediciton set.
                                        combinations)
                
                ## Move prediciton to original coordinates
                prediction <- descale(prediction, mus$curr_var, sigs$curr_var)
                
                ## Store it in the matrix
                predictions[lib_ind, ] <- prediction
                
            } ## Closes for( smp in 1:n_samp )
            
        } ## Closes for (lib_size in lib_sizes)

        ## Save data
        write.table(predictions,
                    file = paste0("runs/", method, "_predictions_", curr_var, ".csv"),
                    quote = FALSE,
                    na = "NA",
                    row.names = FALSE,
                    col.names = FALSE,
                    sep = "," )
                
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
         file = paste0("runs/", method, "_parameters.Rdata") )
    
}

system("rm -f run/*")

save_predictions(file = "originals/three_species.csv",
                 variables = c( "y" ),
                 E = 3, ## Embedding dimension of the system.
                 n_lags = 3, ## 0, -1,..., -(n_lags-1)
                 n_samp = 5, ## Number of random libraries, should be in the hundreds
                 lib = c(501,2001),  ## Library set.
                 pred = c(2501,3000), ## Prediciton set.
                 lib_sizes = (2:4)*10,
             	 method = "we" ## Library sizes
            )

