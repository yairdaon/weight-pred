#!/usr/bin/Rscript
library(rEDM)
source( "helpers/helper.r" )

get_blooms   <- function( df, threshold, ran )
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
    bloom_df <- bloom_df[ range_indices( bloom_df, ran ), ]
    
    return( bloom_df )
}

same_prediction <- function(lib_df,
                            pred_df,
                            method,
                            theta)
{
    ## Stack both data frames 
    stacked <- rbind( lib_df, pred_df )

    ## Probably don't need to sort etc. but just in case...
    stacked <- stacked[ order(stacked$serial_day), ]
    row.names( stacked ) <- 1:nrow( stacked )
    stopifnot( all(
        diff( stacked$serial_day ) > 0
    ) )
    
    ## Load the lib_start, lib_end, pred_start, pred_end variables
    load("data/oos_serial_days.Rdata" )

    lib  <- get_range( stacked,  lib_serial_days )
    pred <- get_range( stacked, pred_serial_days )
    
    print( paste0( "Library = ( ",  lib[1], ", ",  lib[2], " ). Library size = ",  lib[2] -  lib[1] + 1  ) )
    print( paste0( "Test    = ( ", pred[1], ", ", pred[2], " ). Test size = ",    pred[2] - pred[1] + 1  ) )
    print( paste0( "Total number of rows in data frame = ", nrow(stacked) ) )
    
    output <- block_lnlp(stacked,
                         lib = lib,
                         pred = pred,
                         method = method,
                         tp = 0,             
                         columns = c("chl", "silicate_1wk", "AvgDens_1wk", "silicate"),
                         target_column = "chl_p1wk",
                         theta = theta,
                         num_neighbors = -1,
                         stats_only = TRUE,
                         first_column_time = TRUE )
    
    print( paste0( "Rho using best 4D model : ", output$rho ) )
}
        
## Load the data: df, chl_threshold and norm_threshold.
load( "data/processed_block.Rdata" )

## Load the lib_start, lib_end, pred_start, pred_end variables
load("data/oos_serial_days.Rdata" )

full_lib  <- df[ range_indices( df,  lib_serial_days), ]
full_pred <- df[ range_indices( df, pred_serial_days), ]

bloom_lib  <- get_blooms( df, norm_threshold,  lib_serial_days )
bloom_pred <- get_blooms( df, norm_threshold, pred_serial_days )

method <- "s-map"
theta  <- 8
same_prediction(full_lib,
                full_pred,
                method,
                theta)

same_prediction(full_lib,
                bloom_pred,
                method,
                theta )

same_prediction(bloom_lib,
                full_pred,
                method,
                theta )

same_prediction(bloom_lib,
                bloom_pred,
                method,
                theta)

