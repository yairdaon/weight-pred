#!/usr/bin/Rscript
library(rEDM)
source( "plotting.r" )
source( "helper.r" )

lib_start  <- date2serial( "1983-01-01" )
lib_end    <- date2serial( "2010-12-28" )
pred_start <- date2serial( "2011-01-01" )
pred_end   <- date2serial( "2012-03-30" )

###################################
## Prepare the data ###############
###################################

## Read the normalized block
orig_df <- read.csv( "data/block.csv", header = TRUE )

## Reorder the data frame so that , just in case
df <- orig_df[ order(orig_df$serial_day), ]
row.names( df ) <- 1:nrow( df )

## Make sure serial_day increase in time
stopifnot( all(
    diff( df$serial_day ) > 0
) )

args <- commandArgs(trailingOnly = TRUE)
if( grepl("cv", args, ignore.case = TRUE ) )
{
    print( "Leave-one-out cross-validation" )
    ## In leave-one-out-cross-validation we force the library and the
    ## prediction sets to be the same. This tells the code to use
    ## LOOCV.
    lib  <- get_range( df,  lib_start, pred_end )
    pred <- get_range( df,  lib_start, pred_end )
    xtension <- "LOOCV"
} else {
    print( "Out-of-sample results" )
    ## If we do not want to do LOOCV we show out-of-sample skill.
    lib  <- get_range( df,  lib_start,  lib_end )
    pred <- get_range( df, pred_start, pred_end )  
    xtension <- "OoS"
}
ALL <- grepl("all", args, ignore.case = TRUE ) 

n        <- pred[2] - pred[1] + 1

## Start at 4, so we skip serial_day, chl_p1wk and chl
other_vars <- names(df)[4:length(names(df))]
n_vars <- length( other_vars )
combinations <- combn( n_vars, 3 )
num_comb     <- ncol(combinations)

pred_table <- matrix( NA,  nrow = num_comb, ncol = n )
var_table  <- matrix( Inf, nrow = num_comb, ncol = n )

for( i in 1:num_comb )
{
    comb <- combinations[,i]
    ## cols <- c( "chl", combinations[,i] )
    cols <- c( 3, comb + 3 )
      
    output <- block_lnlp(df,
                         lib = lib,
                         pred = pred,
                         method = "s-map",
                         tp = 0, # Time shift is built into the target
                         columns = cols, ##c( 3, comb + 3 ),
                         target_column = 2, ##"chl_p1wk",
                         theta = 8,
                         num_neighbors = -1, ## Take all neighbours
                         first_column_time = FALSE,
                         short_output = TRUE,
                         stats_only = FALSE )

    ## Get model output: prediction times, predictions and variances:
    time <- output[[1]]$model_output$time - pred[1] + 1
    pr   <- output[[1]]$model_output$pred
    vars <- output[[1]]$model_output$pred_var
    
    ## Sanity checks #####################
    len_time <- length(time)
    len_vars <- length(vars)
    len_pred <- length(pr)
    
    stopifnot( len_time == len_vars )
    stopifnot( len_time == len_pred )

    no_vars <- is.na(vars) | is.nan(vars)
    no_pr   <- is.na(pr) 
    if( any(no_vars!=no_pr) )
        stop( "No prediction indices do not align with no variance indices!" )
    ind <- no_pr
    
    ## Save ALL data.
    pr[ ind ] <- NA
    vars[ ind ] <- Inf
    pred_table[ i, time ] <- pr
    var_table [ i, time ] <- vars

}

save(pred_table,
     var_table,
     combinations,
     other_vars,
     n_vars,
     file = paste0("data/", tolower(xtension), "_full_run.Rdata" )
     )
