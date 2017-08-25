#!/usr/bin/Rscript
library(rEDM)
source( "helper.r" )

## Use var to *infer* Chlorophyll data. Should give (?) rho close to
## zero.
var_xmap_chl <- function(var,
                         tp,
                         E=4 )
{
    if( tp > 0 )
        stop( paste0(
            "Do not expect to *predict* Chl from one exogenous variable you can only expect to *infer* it. Take Tp <= 0."
        ) )
    df <- read.csv( "data/raw.csv", header = TRUE )

    ## Extract Chlorophyll time series and insert in the lagged data frame.
    lagged <- data.frame( chl = df$chl )


    ## The target's time series
    var_ts <- df[ , var ]
    
    ## Lag Variable data E more times, for a total of E times.
    for( e in 0:(E-1) )
            lagged[paste0(var,"_", e )] <- half_week_steps( var_ts, tp )

    output <- block_lnlp(lagged,
                         method  = "simplex",
                         columns = 2:(E+1), ## 
                         tp = 0, ## No time step, it is included already.
                         target_column = 1,
                         stats_only = TRUE,
                         first_column_time = FALSE, ## No time variable
                         silent = FALSE)

    
    return( output )
    
}



## Use Chlorophyll data to *infer* var, exogenous variables.
chl_xmap_var <- function(var,
                         tp, ## <= 0 !!!
                         E=4 ## From univariate analysis
                         )
{
    ## Preparations
    if( tp > 0 )
        stop( paste0("Do not expect to *predict* ", var, ", you can only expect to *infer* it. Take Tp <= 0." ) )
    df <- read.csv( "data/univariate_data.csv", header = TRUE )
    
    ## Extract Chlorophyll time series.
    chl_ts <- df$chl
    
    ## The lagged data frame.
    lagged <- data.frame( chl_0 = chl_ts )

    ## Lag Chlorophyll E-1 more times, for a total of E times.
    for( e in 1:(E-1) )        
        lagged[paste0("chl_", e )] <- shift( chl_ts, -e )

    ## The target's time series
    var_ts <- df[ , var ]

    ## Lag the inferred variable. 
    lagged$target <- shift( var_ts, tp )

    ## print( lagged[ !is.na(lagged$target) , ] )
    output <- block_lnlp(lagged,
                         method  = "simplex",
                         columns = 1:E, ## Chlorophyl cols on the left.
                         tp = 0, ## No time step, it is included already.
                         target_column = E+1,
                         stats_only = TRUE,
                         first_column_time = FALSE, ## No time variable
                         silent = FALSE)

    
    return( output )
           
}

