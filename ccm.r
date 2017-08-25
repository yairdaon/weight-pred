#!/usr/bin/Rscript
library(rEDM)
source( "helper.r" )

## Use var to *infer* Chlorophyll data. Should (?) give rho close to
## zero.
var_xmap_chl <- function(var,
                         tp,
                         libsizes,
                         n_samples=100,
                         E=4 )
{
    if( tp > 0 )
        stop( paste0(
            "Do not expect to *predict* Chl from one exogenous variable you can only expect to *infer* it. Take Tp <= 0."
        ) )
    df <- read.csv( "data/raw.csv", header = TRUE )
    
    ## Extract Chlorophyll time series and insert in the lagged data frame.
    lagged <- data.frame( target = half_week_steps( df$chl, tp ) )
    colnames(lagged) <- paste0("chl_", -tp)
    
    ## The target's time series
    var_ts <- df[ , var ]
    
    ## Lag Variable data E more times, for a total of E times.
    for( e in 0:(E-1) )
            lagged[paste0(var,"_", e )] <- half_week_steps( var_ts, -e )
        
    return(
        avg_rho( lagged, n_samples, libsizes )
    )
}



## Use Chlorophyll data to *infer* var, exogenous variables.
chl_xmap_var <- function(var,
                         tp, 
                         libsizes,
                         n_samples=100,
                         E=4 ## From univariate analysis
                         )
{
    ## Preparations
    if( tp > 0 )
        stop( paste0("Do not expect to *predict* ", var, ", you can only expect to *infer* it. Take Tp <= 0." ) )

    df <- read.csv( "data/univariate_data.csv", header = TRUE )
    
    ## The target's time series
    var_ts <- df[ , var ]

    ## Lag the inferred variable. 
    lagged <- data.frame( target = shift( var_ts, tp ) )
    colnames(lagged) <- paste0( var,"_", -tp)
    
    ## Extract Chlorophyll time series.
    chl_ts <- df$chl
    
    ## Lag Chlorophyll E-1 more times, for a total of E times.
    for( e in 0:(E-1) )        
        lagged[paste0("chl_", e )] <- shift( chl_ts, -e )

    return(
        avg_rho( lagged, n_samples, libsizes )
    )               
}


avg_rho <- function(df,
                    n_samples,
                    libsizes)
{
    
    rhos <- numeric( length(libsizes) )
    j <- 1
    
    for( libsize in libsizes )
    {

        for( i in 1:n_samples )
        {
            big    <- sample( nrow(df), libsize + 200, replace = FALSE )
            lib    <- big[ 1 : libsize]
            pred   <- big[ (libsize+1): length(big) ]
            output <- block_lnlp(df,
                                 method  = "simplex",
                                 columns = 2:length( names(df) ),
                                 tp = 0, ## No time step, it is included already.
                                 target_column = 1,
                                 lib = lib,
                                 pred = pred,
                                 stats_only = TRUE,
                                 first_column_time = FALSE, ## No time variable
                                 silent = FALSE)
        
        rhos[j] <- rhos[j] + output$rho
        }
        
        j <- j+1
       
    }
    
    return( rhos / n_samples)    
}

