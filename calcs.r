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

## ## Library and prediction data frames. Not strictly necessary, just
## ## making sure.
## lib_df  <- orig_df[ date_range_indices( orig_df,  lib_start,  lib_end ), ]
## pred_df <- orig_df[ date_range_indices( orig_df, pred_start, pred_end ), ]

## ## Stack both data frames 
## df <- rbind( lib_df, pred_df )
## row.names( df ) <- 1:nrow( df )

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

print( paste0( "Library = ( ",  lib[1], ", ",  lib[2], " ). Library size = ",  lib[2] -  lib[1] + 1  ) )
print( paste0( "Test    = ( ", pred[1], ", ", pred[2], " ). Test size = ",    pred[2] - pred[1] + 1  ) )
print( paste0( "Total number of rows in data frame = ", nrow(df) ) )

all_obs  <- df$chl_p1wk[pred[1]:pred[2]]
all_time <- pred[1]:pred[2]
n        <- length( all_time )
stopifnot( length(all_obs) == n )

## Start at 4, so we skip serial_day, chl_p1wk and chl
other_vars <- names(df)[4:length(names(df))]
n_vars <- length( other_vars )

combinations <- combn( n_vars, 3 )
num_comb     <- ncol(combinations)

exp_vars  <- matrix( 0,   n_vars,               n )       
precision <- matrix( 0,   n_vars,               n )       

sum_exp  <- rep( 0 , n )
sum_var  <- rep( 0 , n )
sum_prec <- rep( 0 , n )
tot_pred <- rep( 0 , n )
cum_pred <- rep( 0 , n )


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
    
    
    ## Keep count of how many predictions were not NA at every single
    ## time.
    tot_pred[time] <- tot_pred[time] + !ind
       
    w <- exp(-vars)
    w[ ind ] <- 0
    sum_exp[time] <- sum_exp[time] + w
    
    u <- 1 / vars
    u[ ind ] <- 0
    sum_prec[time] <- sum_prec[time] + u

    vars[ ind ] <- 0
    sum_var[time] <- sum_var[time] + vars 


    ## Weigh new predictions compare to previous ones
    pr[ ind ] <- 0
    cum_pred[time] <- cum_pred[time] + w * pr  
    
    
}

sum_exp[ sum_exp == 0 ] <- NA
sum_var[ sum_var == 0 ] <- NA
sum_prec[ sum_prec == 0 ] <- NA
predictions <- cum_pred / sum_exp 

filename <- paste0( "data/",xtension, "data.Rdata" ) 
df <- data.frame(serial_day = all_time,
                 chl = all_obs,
                 exp_var = sum_exp / tot_pred,
                 variance = sum_var / tot_pred,
                 precision = sum_prec / tot_pred,
                 predictions = predictions )
save(df,file=filename)


print( paste0( "Skill using weighted predictors: ",  cor(predictions,  all_obs, use = "complete.obs") ) )
