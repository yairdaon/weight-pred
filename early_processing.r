library( zoo )

## Averages ts1 and ts2 (if supplied) over the past week
past_week <- function( ts1, ts2 = NA )
{
    mean_1 <- rollapply(ts1,
                        width = 7,
                        align = 'right',
                        FUN = mean,
                        fill = NA,
                        na.rm = TRUE )
 
    if( is.vector( ts2 ) && length(ts1) == length(ts2) ) {
        mean_2 <- rollapply(ts2,
                            width = 7,
                            align = 'right',
                            FUN = mean,
                            fill = NA,
                            na.rm = TRUE )
        return( ( mean_1 + mean_2 ) / 2 )
    }
    
    return( mean_1 )
}

transform <- function( ts, func = "" )
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
    else if (func == "sqrt")
        ts <- sqrt(ts)
    

    ## Mad is media absolte deviation, scaled so that it is a
    ## consistent estimator of the standard deviation
    ts <- ( ts - median(ts, na.rm = TRUE) ) / mad ( ts, na.rm = TRUE )

    ## If you really really want..
    ## ts <- ( ts - mean(ts, na.rm = TRUE) ) / sd( ts, na.rm = TRUE )
    
    return( ts )
}


## Takes (1, 3, 4, 8, 4, 3, 9) to (NA, 1, 3, 4, 8, 4, 3) 
lag_ts <- function( ts, t = 1 ) {
    if ( t > 0 ) 
        return( c( NA * 1:t, ts[1:(length(ts)-t)] ) )
    else if (t < 0)
        return( c( ts[(1-t):length(ts)], NA * 1:(-t)  ) )
    return( ts )
}

## Returns a time series aligned with the original data frame that has
## a value 6,7 or 8 days ahead (if tp == 2, whichever exists and
## averaged in the unlikely case 6 and 8 exist) or "3.5" days ahead (if
## tp == 1, using either 3 or 4 days ahead, whichever exists).
look_ahead <- function(df,
                       predictee,
                       tp = 1)
{
    if( tp == 1 ) {
        width <- 2
        jump  <- 3 
    } else if( tp == 2 ){
        width <- 3
        jump  <- 6
    } 

    ts  <- df[ ,predictee]
    ts  <- lag_ts( ts, -jump )
    look_ahead <- rollapply(ts,
                            width = width,
                            align = 'left',
                            FUN = mean,
                            fill = NA,
                            na.rm = TRUE )

    return( look_ahead )
}
