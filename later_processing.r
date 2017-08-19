library( zoo )

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

## Take data frame and return only the rows that correspond to the
## serial_days in Hao's dataset.
restrict_day <- function(df, file_name="data/hao_days.txt" )
{
    restriction   <- scan(file_name) 
    restricted_df <- subset(df, serial_day %in% restriction )
    return( restricted_df )
}
