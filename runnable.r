#!/usr/bin/Rscript
library(rEDM)
library(zoo)
source( "helpers/helper.r" )
source( "helpers/mve.r" )

## If you want to download the file for some reason... don't!
## filename = paste0("http://science.sciencemag.org/highwire/filestream/683325/",
##                   "field_highwire_adjunct_files/1/aag0863_SupportingFile_Suppl._Excel_seq1_v2.xlsx")

save_predictions <- function(filename = stop("File name must be provided!"),
                             variables = NULL, 
                             E = 3, ## Embedding dimension of the system.
                             n_lags = E, ## 0,-1, ..., -n_lags
                             n_samp = 150, ## Number of random libraries, should be in the hundreds
                             lib = c(501:2001),  ## Library set.
                             pred = c(2501,3000), ## Prediciton set.
                             lib_sizes = (1:15)*10,     ## Library changes in size and is also random
                             method = "mve"
                             )
{
    ## Save the parameters so we know what parameters we were using
    save(filename,
         variables,
         E,
         n_lags,
         n_samp,
         lib,
         pred,
         lib_sizes,         
         file = paste0("model/runs/", method, "_parameters.Rdata") )
    
    ## Load data 
    raw_df <- read.csv(filename,
                       header = TRUE,
                       sep = "," )

    if( is.null( variables ) )
        variables <- names( raw_df )
    
    pred_func <- mve
    if( method == "uwe" )
        pred_func <- uwe
    
    ## Rescale the data frame and keep track of the normalizing
    ## factors.
    mus  <- get_mus(raw_df)
    sigs <- get_sigs(raw_df)
    df   <- data.frame(scale(raw_df))    
    
    ## Create lags for every variable and make y_p1
    df <- lag_every_variable(df, n_lags)
        
    ## Create the list of combinations
    combinations <- make_combinations(ncol(raw_df), ## TOTAL number of variables
                                      n_lags,
                                      E)
    ## Preallocate
    rhos <- numeric( n_samp )
    
    ## Do the analysis for every variable. 
    for( curr_var in variables )
    {
        empty_file(lib_sizes,
                   filename = paste0("model/runs/", method, "_", curr_var, ".csv") )
                    
        for (lib_size in lib_sizes)
        {
            
            ## For every random library starting point
            for( smp in 1:n_samp )
            {
                if( smp %% 5 == 0 )
                    print( paste0("Var ", curr_var, ", sample ", smp, "/", n_samp, ", lib size ", lib_size ) )
        
                ## Choose random lib
                rand_lib <- random_lib(lib, lib_size)
        
                ## Find the MVE prediction
                prediction <- pred_func(df, ## lagged and scaled.
                                        curr_var,
                                        rand_lib, ## Library set.
                                        pred, ## Prediciton set.
                                        combinations)
                
                ## Move prediciton to original coordinates
                rhos[smp] <- attr( prediction, "rho" )
                
            } ## Closes for( smp in 1:n_samp )
            
            vec <- matrix( c( lib_size, mean(rhos), quantile(rhos, probs = c(0.25,0.5,0.75)) ), nrow = 1)
            write.table(vec,
                        file = paste0("model/runs/", method, "_", curr_var, ".csv"),
                        sep = ",",
                        append = TRUE, 
                        quote = FALSE,
                        col.names = FALSE,
                        row.names = FALSE)                    

        } ## Closes for (lib_size in lib_sizes)
                
    } ## closes for( curr_var in variables )
    
}



## Clean shit up
system("rm -f model/runs/*")

args <- commandArgs( trailingOnly = TRUE )
if( length( args ) > 0 )
{
    save_predictions(file = "model/originals/three_species.csv",
                     variables = c( "y" ),
                     E = 2, ## Embedding dimension of the system.
                     n_lags = 2, ## 0, -1,..., -(n_lags-1)
                     n_samp = 3, ## Number of random libraries, should be in the hundreds
                     lib = c(501,2001),  ## Library set.
                     pred = c(2501,2505), ## Prediciton set.
                     lib_sizes = (2:4)*20,
                     method = "uwe"
                     )
    
    save_predictions(file = "model/originals/three_species.csv",
                     variables = c( "y" ),
                     E = 2, ## Embedding dimension of the system.
                     n_lags = 2, ## 0, -1,..., -(n_lags-1)
                     n_samp = 3, ## Number of random libraries, should be in the hundreds
                     lib = c(501,2001),  ## Library set.
                     pred = c(2501,2505), ## Prediciton set.
                     lib_sizes = (2:4)*20, ## Library sizes
                     method = "mve"
                     )
} else {
    
    save_predictions(file = "model/originals/three_species.csv",
                     variables = c( "y" ),
                     n_samp = 150,
                     lib_sizes = (1:20)*5,
                     method = "mve" 
                     )
    
    save_predictions(file = "model/originals/three_species.csv",
                     variables = c( "y" ),
                     n_samp = 150,
                     lib_sizes = (1:20)*5,
                     method = "uwe" 
                     )
}
