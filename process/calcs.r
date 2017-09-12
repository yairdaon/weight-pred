#!/usr/bin/Rscript
library(rEDM)
source( "helpers/helper.r" )

###################################
## Prepare the data ###############
###################################

## Read the data. Contains df (the data frame), as well as
## chl_threshold and norm_threshold.
load( "data/processed_block.Rdata" )

args <- commandArgs(trailingOnly = TRUE)
if( (length(args) == 0)  ) {
    xtension <- "Test"
} else if (grepl("oos", args, ignore.case = TRUE )) {
    xtension <- "OoS"
} else {
    xtension <- "LOOCV"
}
load(paste0("data/", tolower(xtension), "_serial_days.Rdata"))



lib  <- get_range( df,  lib_serial_days )
pred <- get_range( df, pred_serial_days )  
n    <- pred[2] - pred[1] + 1

## Start at 4, so we skip serial_day, chl_p1wk and chl
other_vars <- names(df)[4:length(names(df))]
n_vars <- length( other_vars )
combinations <- combn( n_vars, 3 )

pred_table <- matrix( NA,  nrow = ncol(combinations), ncol = n )
var_table  <- matrix( Inf, nrow = ncol(combinations), ncol = n )

for( i in 1:ncol(combinations) )
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

chl            <- df$chl       [pred[1]:pred[2]]
chl_p1wk       <- df$chl_p1wk  [pred[1]:pred[2]]
serial_day     <- df$serial_day[pred[1]:pred[2]]
norm_threshold <- quantile(df$chl, 0.95, na.rm = TRUE ) 
save(pred_table,
     var_table,
     combinations,
     other_vars,
     n_vars,
     chl,
     chl_p1wk,
     serial_day,
     norm_threshold,
     file = paste0("data/", tolower(xtension), "_full_run.Rdata" )
     )
