#!/usr/bin/Rscript
library(rEDM)
source( "helper.r" )

date2serial <- function( date_string )
    return(
        as.numeric(as.POSIXlt(as.Date(date_string)))/86400 + 719529
    )

## Gives indices of the rows of df that have dates between dt1 and dt2
date_range_indices <- function( df, start, finish )
    return(
    (df$serial_day >= start) & (df$serial_day <= finish)
    )

get_blooms   <- function( df, threshold, start_date, end_date)
{
    xtend       <- function(n) return( (-7:7) + n )
    
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
    bloom_df <- bloom_df[ date_range_indices( bloom_df, start_date, end_date ), ]
    
    return( bloom_df )
}

same_prediction <- function(lib_df,
                            pred_df,
                            lib_start,
                            lib_end,
                            pred_start,
                            pred_end,
                            method,
                            theta)
{
    ## Stack both data frames 
    stacked <- rbind( lib_df, pred_df )

    stacked <- stacked[ order(stacked$serial_day), ]
    row.names( stacked ) <- 1:nrow( stacked )
    stopifnot( all(
        diff( stacked$serial_day ) > 0
    ) )
    
    lib_range <- which( date_range_indices( stacked, lib_start, lib_end ) )
    lib <- c( min(lib_range), max(lib_range) )

    pred_range <- which( date_range_indices( stacked, pred_start, pred_end ) )
    pred <- c( min(pred_range), max(pred_range) )
    
    print( paste0( "Library = ( ", lib[1],  ", ", lib[2],  " )"   ) )
    print( paste0( "Test    = ( ", pred[1], ", ", pred[2], " )"   ) )
    print( paste0( "Total number of rows in data frame = ", nrow(df) ) )
    
    output <- block_lnlp(df,
                         lib = lib,
                         pred = pred,
                         method = method,
                         tp = 0,                    
                         columns = c("chl", "silicate_m1wk", "AvgDens_1wk", "silicate"),
                         target_column = "chl_p1wk",
                         theta = theta,
                         stats_only = TRUE )
    
    print( paste0( "OOS rho using best 4D model : ", output$rho ) )
}
        
normalize <- function( ts, m, s )
    return( (ts - m)/s )

    
## Prepare the data.
load( "originals/chl_block_full.Rdata" )

chl_mean <- mean( orig_block$chlA,        na.rm = TRUE )
chl_sd   <- sd(   orig_block$chlA,        na.rm = TRUE )
sil_mean <- mean( orig_block$silicate,    na.rm = TRUE )
sil_sd   <- sd(   orig_block$silicate,    na.rm = TRUE )
den_mean <- mean( orig_block$AvgDens_1wk, na.rm = TRUE )
den_sd   <- sd(   orig_block$AvgDens_1wk, na.rm = TRUE )

df <- data.frame(
    serial_day    = orig_block$serial_day,
    chl_p1wk      = normalize( orig_block$chlA_.1wk,    chl_mean, chl_sd ),
    chl           = normalize( orig_block$chlA,         chl_mean, chl_sd ),
    silicate      = normalize( orig_block$silicate,     sil_mean, sil_sd ),
    silicate_m1wk = normalize( orig_block$silicate_1wk, sil_mean, sil_sd ),
    AvgDens_1wk   = normalize( orig_block$AvgDens_1wk,  den_mean, den_sd )
)
row.names( df ) <- 1:nrow(df)

threshold  <- quantile( df$chl, 0.95, na.rm = TRUE )

lib_start  <- date2serial( "1983-01-01" )
lib_end    <- date2serial( "2010-12-28" )
pred_start <- date2serial( "2011-01-01" )
pred_end   <- date2serial( "2012-03-30" )


full_lib  <- df[ date_range_indices( df,  lib_start,  lib_end ), ]
full_pred <- df[ date_range_indices( df, pred_start, pred_end ), ]

bloom_lib  <- get_blooms( df, threshold,  lib_start,  lib_end )
bloom_pred <- get_blooms( df, threshold, pred_start, pred_end )

method <- "s-map" ##"simplex"
theta <- 8
same_prediction(full_lib,
                full_pred,
                lib_start,
                lib_end,
                pred_start,
                pred_end,
                method,
                theta)

same_prediction(full_lib,
                bloom_pred,
                lib_start,
                lib_end,
                pred_start,
                pred_end,
                method,
                theta )

same_prediction(bloom_lib,
                full_pred,
                lib_start,
                lib_end,
                pred_start,
                pred_end,
                method,
                theta )

same_prediction(bloom_lib,
                bloom_pred,
                lib_start,
                lib_end,
                pred_start,
                pred_end,
                method,
                theta)


