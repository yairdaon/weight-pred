#!/usr/bin/Rscript
library(rEDM)
library(zoo)
source( "helpers/helper.r" )
source( "helpers/mve.r" )

## If you want to download the file for some reason... don't!
## filename = paste0("http://science.sciencemag.org/highwire/filestream/683325/",
##                   "field_highwire_adjunct_files/1/aag0863_SupportingFile_Suppl._Excel_seq1_v2.xlsx")

save_predictions <- function(filename = stop("File name must be provided!"),
                             target_column = "y", 
                             E = 3, ## Embedding dimension of the system.
                             max_lag = E, ## 0,-1, ..., -max_lag
                             n_samp = 100, ## Number of random libraries, should be in the hundreds
                             lib = c(501:2000),  ## Library set.
                             pred = c(2500,2999), ## Prediciton set.
                             lib_sizes = (1:25)*5,    ## Library changes in size and is also random
                             method = "mve",
                             num_neighbors = E+1 )
{
    print( paste0( "Method: ", method, ", target: ", target_column) )
    
    ## Load data 
    raw_df <- read.csv(filename,
                       header = TRUE,
                       sep = "," )
    
    ## ## Rescale the data frame and keep track of the normalizing
    ## ## factors.
    mus  <- get_mus(raw_df)
    sigs <- get_sigs(raw_df)
    df   <- data.frame(scale(raw_df))    

    filename <- paste0("model/runs/", method, "_", target_column, "_", num_neighbors,".csv")                
    empty_file( filename = filename )
    
    pred_func <- mve
    if( method == "uwe" )
        pred_func <- uwe
    if( method == "hao" ) 
        pred_func <- multiview

    ## Preallocate
    rhos <- numeric( n_samp )
        
    for (lib_size in lib_sizes)
    {
        
        ## For every random library starting point
        for( smp in 1:n_samp )
        {
            if( smp %% 5 == 0 )
                print( paste0("Sample ", smp, "/", n_samp, ", lib size ", lib_size ) )
            
            ## Find the MVE prediction
            output <- pred_func(df,
                                lib = random_lib(lib, lib_size),
                                pred = pred,
                                ## norm_type = c("L2 norm", "L1 norm", "P norm"),
                                ## P = 0.5, 
                                E = E,
                                tau = 1,
                                tp = 1,
                                max_lag = max_lag,
                                num_neighbors = num_neighbors,
                                k = "sqrt",
                                na.rm = FALSE, 
                                target_column = target_column, 
                                stats_only = TRUE,
                                first_column_time = FALSE, 
                                exclusion_radius = NULL,
                                silent = FALSE)
            rhos[smp] <- output$rho
            
        } ## Closes for( smp in 1:n_samp )
                
        ## Order agrees with order set in empty_file
        vec <- matrix( c( lib_size, mean(rhos), quantile(rhos, probs = c(0.25,0.5,0.75)) ), nrow = 1)
        write.table(vec,
                    file = filename, 
                    sep = ",",
                    append = TRUE, 
                    quote = FALSE,
                    col.names = FALSE,
                    row.names = FALSE)                    
        
    } ## Closes for (lib_size in lib_sizes)

    ## Save the parameters so we know what parameters we were using
    ## save(filename,
    ##      E,
    ##      max_lag,
    ##      n_samp,
    ##      lib,
    ##      pred,
    ##      lib_sizes,         
    ##      file = "model/runs/parameters.Rdata" )
    
} ## Closes function save_predictions

args <- commandArgs( trailingOnly = TRUE )
if( length(args) > 2 && args[3] == "test" )
{
    print( "Testing..." )
    n_samp <- 5
    lib_sizes <- c(25,50)
} else {
    n_samp <- 100
    lib_sizes <- (1:10)*10
}




save_predictions(file = "model/originals/three_species.csv",
                 target_column = "y",
                 n_samp = n_samp,
                 method = args[1],
                 num_neighbors = as.numeric(args[2]),
                 lib_sizes = lib_sizes, 
                 )
    
