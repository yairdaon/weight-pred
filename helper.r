library( zoo )

hold <- function(k=0){
    message("Press Return To Continue")
    invisible(readLines("stdin", n=k))
}

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

harry_plotter <- function(x, tit, k=0 )
{
    tit <- paste( "Transformed ", tit, " data" )
    x <- na.omit( x )
    x11()
    hist(x,
         main = tit,
         breaks = 20,
         probability=TRUE,
         xlim = c(-4, 4 ) )
    
    lines( density(x), col="red" )
    z <- seq(-4,4, length = 1000)
    lines(z,
          dnorm(z, mean=0, sd = 1 ),
          col = "blue" )
    if( k == 0 )
        hold()
    else
        hold(k)
    
}

## Takes (1, 3, 4, 8, 4, 3, 9) to (NA, 1, 3, 4, 8, 4, 3) 
lag_ts <- function( ts ) {
    return( c( NA, ts[1:length(ts)-1] ) )
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
