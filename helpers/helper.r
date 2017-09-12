library( zoo )

## Gives indices of the rows of df that have dates between dt1 and dt2
range_indices <- function( df, ran )
    return(
    (df$serial_day >= ran[1]) & (df$serial_day <= ran[2])
    )


get_range <- function( df, ran )
{
    range <- which( range_indices(df,ran) ) 
    return( c( min(range), max(range) ) )
}

transform <- function(ts,
                      func = "",
                      robust = TRUE )
{
    minimum <- min(ts, na.rm = TRUE)
    
    ## If minimum is zero or less ...
    if (  minimum <= 0 && (func == "log" || func == "sqrt") ) {

        ## ... make sure the minimum is zero...
        ts <- ts - minimum
        
        ## ... and make sure the minimum is actually half the second
        ## smallest value (which is guaranteed positive).
        ts[ ts == 0 ]  <- min( ts[ts != 0], na.rm = TRUE )/2
           
    }

    
    ## Transform data according to a given function.
    if( func == "log" )
        ts <- log( ts )
    else if ( func == "sqrt" )
        ts <- sqrt(ts)
    

    if( robust )
        ## Mad is media absolte deviation, scaled so that it is a
        ## consistent estimator of the standard deviation
        ts <- ( ts - median(ts, na.rm = TRUE) ) / mad ( ts, na.rm = TRUE )
    else
        ## If you really really want to use mean and std-dev...
        ts <- ( ts - mean(ts, na.rm = TRUE) ) / sd( ts, na.rm = TRUE )
        
        return( ts )
}

