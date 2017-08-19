library( zoo )

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


## Returns a data frame with lags of time series ts, named var_0,
## var_1 etc. lagged E-1 times.
multi_lag_ts <- function(dat, var, E) 
{  
    df <- data.frame( dat$serial_day, dat[ ,var ])
    colnames( df ) <- c( "serial_day", paste0(var, "_0") )
    
    ## Predictor at time t-E
    for( e in 1:(E-1) )
        df[paste0(var, "_", e )] <- lag_ts( df[,paste0(var, "_", e-1 )] )

    ## print( df[1:10,] )
    return( df )
}

legal_indices <- function(serial_day,
                          numLags,
                          which_obs )
{
    diffVec  <- c( NA, diff(serial_day) )
    indices0 <- TRUE
    indices1 <- TRUE
    indices2 <- TRUE

    ## Calculate indices of valid observations:
    for( E in 1:(numLags-1) )
    {
        ## Observations preceded by half weeks obeservations
        indices0 <- indices0 & ( diffVec < 5 ) & ( diffVec > 2 )
        
        ## Take care of alternating observations, 3-4-3-4-3... and
        ## 4-3-4-3-4...
        if( E %% 2 == 1 ) {
            indices1 <- indices1 & diffVec == 3 
            indices2 <- indices2 & diffVec == 4
        } else {
            indices1 <- indices1 & diffVec == 4 
            indices2 <- indices2 & diffVec == 3
        }
    
        ## Lag the time difference vector
        diffVec <- lag_ts( diffVec )
    }

    ## Exactly alternating observations - Merge the observations above
    indices3 <- indices1 | indices2
    
    ## IF there is a NA it means there wasn't a long enough sequence
    ## of observations and the corresponding observation shoud be
    ## omitted.
    indices0[is.na(indices0)] <- FALSE
    indices1[is.na(indices1)] <- FALSE
    indices2[is.na(indices2)] <- FALSE
    indices3[is.na(indices3)] <- FALSE
    
    ## Keep only the parts of the data frame that are useful,
    ## according to our choice of which observations to keep
    if( which_obs == "34" )
        return( indices1 )
    else if( which_obs == "43" )
        return( indices2 )
    else if( which_obs == "alt" )
        return( indices3 )
    else
        return( indices0 )
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

## Take data frame and return only the rows that correspond to the
## serial_days in Hao's dataset.
restrict_day <- function(df, file_name="data/hao_days.txt" )
{
    restriction   <- scan(file_name) 
    restricted_df <- subset(df, serial_day %in% restriction )
    return( restricted_df )
}
