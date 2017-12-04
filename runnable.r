#!/usr/bin/Rscript
library(rEDM)
library(zoo)
source( "helpers/helper.r" )
source( "helpers/mve.r" )
set.seed( 19 )

## filename = paste0("http://science.sciencemag.org/highwire/filestream/683325/",
##                   "field_highwire_adjunct_files/1/aag0863_SupportingFile_Suppl._Excel_seq1_v2.xlsx")

save_predictions <- function(dirname = stop("Directory name must be provided!"),
                             target_column = 1, 
                             E = 3, ## Embedding dimension of the system.
                             max_lag = E, ## 0,-1, ..., -max_lag
                             n_samp = 100, ## Number of random libraries, should be in the hundreds
                             lib = c(501,2000),  ## Library set.
                             pred = c(2500,2999), ## Prediciton set.
                             lib_sizes = (1:25)*5,    ## Library changes in size and is also random
                             method = "mve",
                             num_neighbors = E+1 )
{
    print( paste0( "Method: ", method, ", target: ", target_column) )
    
    ## Load data 
    raw_df <- read.csv(paste0( dirname, "/original.csv" ),
                       header = TRUE,
                       sep = "," )
    
    ## ## Rescale the data frame and keep track of the normalizing
    ## ## factors.
    mus   <- get_mus(raw_df)
    sigs  <- get_sigs(raw_df)
    df    <- data.frame(scale(raw_df))
    noise <- rnorm( prod(dim(df)), mean=0, sd=sqrt(0.1))  
    df    <- noise + df ## As long as first_column_time == FALSE
    
    filename <- paste0(dirname, "/runs/", method, "_", target_column, "_", num_neighbors,".csv")                
    empty_file( filename = filename )
    
    pred_func <- mve
    if( method == "uwe" )
        pred_func <- uwe
    if( method == "hao" )
    {
        pred_func <- multiview
        target_column <- which( names(raw_df) == target_column )
    }
    ## Preallocate
    rhos <- numeric( n_samp )
        
    for (lib_size in lib_sizes)
    {
        
        ## For every random library starting point
        for( smp in 1:n_samp )
        {            
            rand_lib <- random_lib(lib, lib_size)

            ## Find the MVE prediction
            output <- pred_func(df,
                                lib = rand_lib,
                                pred = pred,
                                norm_type = c("L2 norm", "L1 norm", "P norm"),
                                P = 0.5, 
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
            if( smp %% 5 == 0 )
            {
                avg <- mean(rhos[1:smp])
                err <- sd( rhos[1:smp] )/ sqrt( smp )
                print( paste0("Sample ", smp, "/", n_samp, ", lib size ", lib_size, " mean ", avg, " +/- ", err ) )
                if (err < 0.01*avg )
                    break
            }
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
    
} ## Closes function save_predictions

args <- commandArgs( trailingOnly = TRUE )

if( args[1] == "huisman" ) {
    dirname <-"huisman"
} else if( args[1] == "hp" ) {
    dirname <- "hp"
} else {
    stop( "Unknown model" )
}

if( args[2] == "mve" ) {
    method <- "mve" 
} else if( args[2] == "uwe" ) {
    method <- "uwe"
} else {
    stop( "Unknown method" )
}

target_column <- args[3]
num_neighbors <- as.numeric(args[4] )
  
if( args[length(args)] == "test" ) {
    print( "Testing..." )
    n_samp <- 5
    lib_sizes <- c(25)
} else if( args[3] == "long" ) {
    n_samp <- 50
    lib_sizes <- c(25,50)
} else {
    n_samp <- 150
    lib_sizes <- (1:10)*10
}

save_predictions(dirname = dirname,
                 target_column = target_column,
                 n_samp = n_samp,
                 method = method,
                 num_neighbors = num_neighbors,
                 lib_sizes = lib_sizes, 
                 )
    
