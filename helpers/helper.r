library( zoo )

## Craeate lags for every variable
lag_every_variable <- function(df, n_lags)
{
    for( var in names(df) )
    {
        ts <- zoo( df[ , var ] )
        for( k in 1:(n_lags-1) )
            df[ paste0(var, "_", k) ] <- c(rep( NA, k ),
                                           lag( ts, -k )
                                           )
    }
    return( df )
}


## Guess what this function does...
normalize_df <- function(df)
{
    vars <- names(df) 
    df <- scale(df) ## This makes the df a matrix...
    df <- data.frame( df ) ## ...so make it a df...
    names(df) <- vars ##...and restore the names.
    return( df )
}


bloom_or_not_days <- function(df,
                              threshold,
                              ran,
                              bloom = TRUE)
{
    ## Used later...
    xtend       <- function(n) return( (-14:14) + n )

    ## Get serial days and chlorophyll levels, restricted according to
    ## ran.
    ind <- range_indices(df,ran)
    serial_day <- df$serial_day[ ind ]
    chl        <- df$chl       [ ind ]
    
    ## Find the days above threshold.
    bloom_days <- serial_day[ chl > threshold ]

    ## Extend the bloom days two weeks in both directions
    bloom_days <- unique( as.vector( sapply(bloom_days,xtend) ) )

    ## Keep only the days that are in the data frame
    bloom_days <- bloom_days[ bloom_days %in% serial_day ]

    ## All the rest, by definition.
    no_bloom_days <- serial_day[ !(serial_day %in% bloom_days) ] 
    
    if ( length(bloom_days) + length(no_bloom_days) != length(serial_day) )
        stop( paste0( length(bloom_days), " + ", length(no_bloom_days), " != ", length(serial_day) ) ) 
             
    if( bloom )
        return( bloom_days )
    else
        return( no_bloom_days )
}


get_blooms <- function(df,
                       threshold,
                       ran)
{

    bloom_days <- bloom_or_not_days(df,
                                    threshold,
                                    ran,
                                    TRUE)

    
    ## Restrict the data frame to the above mentioned bloom days
    bloom_df <- df[ df$serial_day %in% bloom_days, ]
    rownames( bloom_df ) <- 1:nrow(bloom_df)
    
    ## Restrict to the desired date range
    bloom_df <- bloom_df[ range_indices( bloom_df, ran ), ]
    
    return( bloom_df )
}

get_no_blooms <- function(df,
                          threshold,
                          ran)
{
    no_bloom_days <- bloom_or_not_days(df,
                                       threshold,
                                       ran,
                                       FALSE)
    
    ## Restrict the data frame to the non-bloom days
    no_bloom_df <- df[ df$serial_day %in% no_bloom_days, ]
    rownames( no_bloom_df ) <- 1:nrow(no_bloom_df)

    return( no_bloom_df )
}

date2serial <- function( date_string )
    return(
        as.numeric(as.POSIXlt(as.Date(date_string)))/86400 + 719529
    )

serial2date <- function( serial_day )
    return(
        as.Date(serial_day - 719529, origin = "1970-01-01") 
    )

## ## Tests, if u wanna use them
## ex_date <- "2008-02-29"
## print( ex_date )
## ex_day <- date2serial( ex_date ) 
## print( ex_day )
## rec_date <- serial2date( ex_day )
## print( rec_date )
## rec_day <- date2serial( rec_date )
## print( rec_day )




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

