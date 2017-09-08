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
    print( "Leave-one-out-cross-validation" )
    ## In leave-one-out-cross-validation we force the library and the
    ## prediction sets to be the same. This tells the code to use
    ## LOOCV.
    lib  <- get_range( df,  lib_start, pred_end )
    pred <- get_range( df,  lib_start, pred_end )
    xtension <- "LOOCV"
    load("data/loocv_full_run.Rdata")
} else {
    print( "Out-of-sample" )
    ## If we do not want to do LOOCV we show out-of-sample skill.
    lib  <- get_range( df,  lib_start,  lib_end )
    pred <- get_range( df, pred_start, pred_end )  
    xtension <- "OoS"
    load("data/oos_full_run.Rdata")
}

all_obs  <- df$chl_p1wk[pred[1]:pred[2]]
all_time <- pred[1]:pred[2]
n        <- length( all_time )
stopifnot( length(all_obs) == n )


## pred_table[ is.na(pred_table) ] <- 0
## exp_var <- exp(-var_table)
## cum <- colSums( exp_var )
## weight_pred <- colSums( exp_var * pred_table ) / cum
## print( paste0( "Skill using weighted predictors: ",  cor(weight_pred, all_obs, use = "complete.obs") ) )


## mins <- apply( var_table, 2, which.min )
## min_vars_pred <- numeric(n)
## for( i in 1:n )
##     min_vars_pred[i] <- pred_table[ mins[i] , i ]
## print( paste0( "Skill using min-var  predictors: ",  cor(min_vars_pred, all_obs, use = "complete.obs") ) )


probs <- (1:1000) / 1000
exp_weight <- matrix(0,
                      nrow = length(probs),
                      ncol = n )

prec_weight <- matrix(0,
                      nrow = length(probs),
                      ncol = n )


for( i in 1:n )
{
    v <- var_table[ , i]
    Q <- quantile( v, probs )
    j <- 1
    for( q in Q  )
    {
        ind <- which( v < q )
        p <- pred_table[ ind, i]
        exp_w  <- exp(-v[ind])
        prec_w <- 1/v[ind] 
        exp_weight [j,i] <- sum(p*exp_w) / sum(exp_w)
        prec_weight[j,i] <- sum(p*prec_w) / sum(prec_w)
        j <- j + 1
    }  
}

exp_rhos  <- numeric(length(probs))
prec_rhos <- numeric(length(probs))
for( j in 1:length(probs) )
{
    exp_rho  <- cor( exp_weight[j, ], all_obs, use = "complete.obs")
    prec_rho <- cor(prec_weight[j, ], all_obs, use = "complete.obs")
    
    exp_rhos[j]  <- exp_rho
    prec_rhos[j] <- prec_rho
}


png(filename=paste0("plots/exp_weighted_predictions_", xtension, ".png") )
plot(probs,
     exp_rhos,
     type = "l",
     xlab = "Quantile",
     ylab = "Rho",
     main = paste0( xtension, " skill per quantile for exp-weights" )
     )
dev.off()

png(filename=paste0("plots/precision_weighted_predictions_", xtension, ".png") )
plot(probs,
     prec_rhos,
     type = "l",
     xlab = "Quantile",
     ylab = "Rho",
     main = paste0( xtension, " skill per quantile for precision-weights" )
     )
dev.off()
