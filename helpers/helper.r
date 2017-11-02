#!/usr/bin/Rscript
library( zoo )

make_combinations <- function(n_vars,
                              n_lags,
                              E)
{
    ## We know the allowed combinations must have one unlagged
    ## variable. By the way function combn is structured, these will be
    ## the first n_comb (defined below).
    combinations <- combn( n_vars*n_lags, E ) ## Get ALL cobminations.
    ## Only this many are OK: subtract the bad ones
    n_comb <- choose( n_vars*n_lags, E )- choose(n_vars*(n_lags-1), E ) 
    ## Throw away the rest.
    combinations <- combinations[ ,1:n_comb ]

    return(combinations)
}
        
get_mus <- function(df)
{
    mus  <- colMeans(df, na.rm = TRUE)
    mus  <- setNames(as.list(mus), names(df))
    return(mus)
}

get_sigs <- function(df)
{
    sigs <- apply(df, 2, sd, na.rm = TRUE )
    sigs <- setNames(as.list(sigs), names(df))
    return(sigs)
}

random_lib <- function(lib_range, lib_size )
{
    lib_start <- sample(lib_range, 1)
    lib       <- c(lib_start, lib_start + lib_size- 1 )
}

## Get a weighted prediction using the predictors and their
## uncertainty.
weighted_prediction <- function(var_table,
                                pred_table)
{
    ## Sanity check
    stopifnot(
        nrow(var_table) == nrow(pred_table) &&
        ncol(var_table) == ncol(pred_table)
    )

    ## Meomry allocation
    ret <- numeric(ncol(pred_table))

    ## Avoid loops at all costs!!! Or be lazy!!
    for( i in 1:ncol(pred_table) )
    {            
        ## Extract uncertainty of all predictors and find its
        v   <- var_table[ , i]
        ind <- order( v )[1:ceiling(sqrt(nrow(var_table)))]
        
        ## Find who is in, get corresponding predictions and get
        ## corresponding weights, both exponential and precision.
        p <- pred_table[ind, i]
        w <- exp(-v[ind])
            
        ## Weight the predictors according to the weights and store
        ## the weighted prediction in the preallocated matrices.
        ret[i] <- sum(p * w) / sum(w)
    }
    return( ret )
}


## Craeate all lags and push-ahead for every variable. If 
lag_every_variable <- function(df, n_lags)
{
    ## Take every variable ...
    for( var in names(df) )
    {
        ## ... make that variable a time series ...
        ts <- zoo( df[ , var ] )

        ## ... and create its lags.
        for( k in 1:(n_lags-1) )
            df[ paste0(var, "_", k) ] <- c(rep( NA, k ),
                                           lag( ts, -k )
                                           )
    }

    ## Create dataframe to hold the look-ahead time series
    tmp <- data.frame(tmp = numeric(nrow(df)))
    tmp$tmp <- NULL
    
    ## Take every variable ...
    for( var in names(df) )
        ## ... and and create a look-ahead time series.
        tmp[ paste0(var, "_p1") ] <- c( as.vector(lag(ts,1)) , NA )

    ## Merge the dataframes so that the names are sorted as follows:
    ## look-ahead variables, non-lagged variables, first variable lags, second variable lags, etc.
    ## E.g. x_p1, y_p1, z_p1, x, y, z, x_1, x_2, y_1, y_2, z_1, z_2
    df <- data.frame( c(tmp,df) )
    ## print( names(df) ) ## If you don't believe what I wrote above.
    
    return( df )
}
## ## Normalize data frame AND keep track of column means and standard
## ## deviations.
## normalize_df <- function(df)
## {
##     ## Keep track of the names
##     vars <- names(df)

##     ## Hold the means and sds
##     means <- colMeans(df, na.rm = TRUE)
##     sds   <- apply(df, 2, sd, na.rm = TRUE )

##     df <- scale(as.matrix(df))
##     df <- data.frame(df)
##     names(df) <- vars

##     df <- (df - means)/sds
##     attr(df, "means") <- means
##     attr(df, "sds") <- sds
##     return(df)
## }


descale <- function( df, mu = NULL, sig = NULL )
{

    ## Undo the scaling
    if( is.null(sig) ) {
        if( !is.null(attr(df,"scaled:scale")) )
            df <- scale(df, center = FALSE, scale = 1/attr(df, "scaled:scale" ) )
    }
    else
        df <- scale(df, center = FALSE, scale = 1/sig)
    
    ## Undo the centering
    if( is.null(mu) ) {
        if( !is.null(attr(df,"scaled:center"))) 
            df <- scale(df, center = -attr(df, "scaled:center"), scale = FALSE  )
    }
    else 
        df <- scale(df, center = -mu, scale = FALSE )
    
    return(df)
}
## Tests for function rescale
df <- data.frame(x = c(1,4,6,2,4,6,-1),
                 y = c(1,2,3,4,5,6,7),
                 z = c(0,1,0,1,0,1,0))
df1 <- scale( df )
df2 <- descale(df1)
stopifnot(all(
    abs(df2-df) < 1e-14
))

df3 <- descale(df1,
               mu  = attr(df1, "scaled:center"),
               sig = attr(df1, "scaled:scale" )  
               )
stopifnot(all(
    abs(df3-df) < 1e-14
))

true_cor <- function(x, y, mu, sig ) {

    x <- descale(x, mu = mu, sig = sig )
    y <- descale(y, mu = mu, sig = sig )

    return( cor(x, y, use = "complete.obs" ) )
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

