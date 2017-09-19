library( zoo )

get_blooms   <- function( df, threshold, ran )
{
    xtend       <- function(n) return( (-17:17) + n )
    
    ## Find the bloom rows
    bloom_rows <- which( df$chl > threshold )

    ## Using the rows, find bloom days
    bloom_days <- df$serial_day[bloom_rows]

    ## Extend the bloom days one week in both directions
    bloom_days <- unique( as.vector( sapply(bloom_days, xtend ) ) )

    ## Restrict the data frame to the above mentioned bloom days
    bloom_df   <- df[ df$serial_day %in% bloom_days, ]
    rownames( bloom_df ) <- 1:nrow(bloom_df)

    ## Restrict to the desired date range
    bloom_df <- bloom_df[ range_indices( bloom_df, ran ), ]
    
    return( bloom_df )
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

