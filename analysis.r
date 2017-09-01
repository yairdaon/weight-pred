#!/usr/bin/Rscript
library(rEDM)
source( "plotting.r" )
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

get_range <- function( df, start, finish )
{
    range <- which( date_range_indices( df, start, finish ) )
    return( c( min(range), max(range) ) )
}

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

print( paste0( "Library = ( ",  lib[1], ", ",  lib[2], " )"   ) )
print( paste0( "Test    = ( ", pred[1], ", ", pred[2], " )"   ) )
print( paste0( "Total number of rows in data frame = ", nrow(df) ) )

all_pred_obs  <- df$chl_p1wk[pred[1]:pred[2]]
all_pred_time <- pred[1]:pred[2]
n             <- length( all_pred_time )
stopifnot( length(all_pred_obs) == n )

########################################
## Prediction using the best 4D model ##
########################################
best <- block_lnlp(df,
                   lib  = lib,
                   pred = pred,
                   method = "s-map",
                   tp = 0,                    
                   columns = c("chl", "silicate_1wk", "AvgDens_1wk", "silicate"),
                   target_column = "chl_p1wk",
                   theta = 8,
                   short_output = TRUE,
                   stats_only = FALSE )

best_pred       <- rep( NA, n )
time            <- best[[1]]$model_output$time - pred[1] + 1
best_pred[time] <- best[[1]]$model_output$pred
best_rho        <- best[[1]]$stat$rho

jpeg( paste0( "plots/best_pred_", xtension, ".jpg" ) )
plot(all_pred_time,
     all_pred_obs,
     type = "l",
     col  = "black",
     main = "Best model predictions compared to truth")
lines(all_pred_time,
      best_pred,
      col = "red" )

###########################################
## Now consider the averaged predictions ##
###########################################
min_vars   <- rep( Inf, n )
var_pred   <- rep( NA,  n )
cum_pred   <- numeric(n)
cum_weight <- numeric(n)

## Explicitly state the variables allowed to be used for prediction
variables <- c("nitrate",
               "nitrate_1wk",
               "nitrate_2wk",
               "silicate",
               "silicate_1wk",
               "silicate_2wk",
               "nitrite",
               "nitrite_1wk",
               "nitrite_2wk",
               "AvgTemp_1wk",
               "AvgTemp_1wk_1wk",
               "AvgTemp_1wk_2wk",
               "AvgDens_1wk",
               "AvgDens_1wk_1wk",
               "AvgDens_1wk_2wk",
               "U_WIND",
               "U_WIND_1wk",
               "U_WIND_2wk",
               "chl_1wk",
               "chl_2wk" )
               ## "chl", NOT including chl because it is included by default
combinations <- combn( variables, 3 )

for( i in 1:dim(combinations)[2] )
{
    cols <- c( "chl", combinations[,i] )

    output <- block_lnlp(df,
                         lib = lib,
                         pred = pred,
                         method = "s-map",
                         tp = 0, # Time shift is built into the target
                         columns = cols,
                         target_column = "chl_p1wk",
                         theta = 8,
                         short_output = TRUE,
                         stats_only = FALSE )

    time <- output[[1]]$model_output$time - pred[1] + 1

    ## The model's predictions
    pr <- rep( NA, n )
    pr[time] <- output[[1]]$model_output$pred

    ## The prediction's associated variances
    vars <- rep( Inf, n )
    vars[ time ] <- output[[1]]$model_output$pred_var
    vars[ is.nan(vars) ] <- Inf
    vars[ is.na (vars) ] <- Inf
    
    ## Sanity checks #####################
    len_time <- length(time)
    len_pred <- length(output[[1]]$model_output$pred)
    len_vars <- length(output[[1]]$model_output$pred_var)

    stopifnot( len_time == len_pred )
    stopifnot( len_time == len_vars )
    
    infys <- (vars==Inf)
    nas   <- is.na(pr) 
    if( any(infys!=nas) )
        stop( "Infinite vairance and NA predictions do not match!" )

    ## Where we don't have prediction set it to some arbitrary
    ## value. This will not propagate to the final prediction since we
    ## force their corresponding variance to be infinite.
    pr[ is.na(pr) ] <- 0 
    stopifnot( !any( is.na(vars) ) )
    
    ## Choose new predictions if they are less uncertain than previous
    ## ones.
    min_ind                   <- which( vars < min_vars )
    min_vars[min_ind]         <- vars[min_ind]
    var_pred[min_ind]         <- pr[min_ind]

    ## Weigh new predictions compare to previous ones
    weight     <- exp(-2*vars)
    cum_pred   <- cum_pred + weight * pr  
    cum_weight <- cum_weight + weight

}

cum_pred <- cum_pred / cum_weight
two_pred <- ( cum_pred + var_pred ) / 2

stopifnot( length(all_pred_time) == length(best_pred   ) )
stopifnot( length(all_pred_time) == length(cum_pred    ) )
stopifnot( length(all_pred_time) == length(var_pred    ) )


jpeg( paste0( "plots/avg_pred_", xtension, ".jpg" ) )
plot(all_pred_time,
     all_pred_obs,
     type = "l",
     col = "black",
     main = "Averaged model predictions compared to truth")
lines(all_pred_time,
      cum_pred,
      col = "red" )

jpeg( paste0( "plots/minvar_pred_", xtension, ".jpg" ) )
plot(all_pred_time,
     all_pred_obs,
     type = "l",
     col  = "black",
     main = "Minimal variance predictions compared to truth")
lines(all_pred_time,
      var_pred,
      col = "blue" )

jpeg( paste0( "plots/both_pred_", xtension, ".jpg" ) )
plot(all_pred_time,
     all_pred_obs,
     type = "l",
     col  = "black",
     main = "Two models predictions compared to truth")
lines(all_pred_time,
      two_pred,
      col = "purple" )


print( paste0( "Skill using best predictor:      ",  best_rho ) )
print( paste0( "Skill using weighted predictors: ",  cor(cum_pred,  all_pred_obs, use = "complete.obs") ) )
print( paste0( "Skill using min var predictors:  ",  cor(var_pred,  all_pred_obs, use = "complete.obs") ) )
print( paste0( "Skill averaging both previous:   ",  cor(two_pred,  all_pred_obs, use = "complete.obs") ) )

## warnings()
