#!/usr/bin/Rscript
library(rEDM)

args <- commandArgs(trailingOnly = TRUE)
if( length(args) == 0 ) {
    xtension <- "Test"
    probs <- (1:20) / 20
} else if( grepl("oos", args, ignore.case = TRUE ) ) {
    xtension <- "OoS"
    probs <- (1:1000) / 1000
} else {
    xtension <- "LOOCV"
    probs <- (1:1000) / 1000
}

load(paste0("data/", tolower(xtension), "_full_run.Rdata"))
dims <- dim( pred_table )
n <- dims[2]

## Preallocate memory, for speed
exp_weight  <- matrix(0, nrow = length(probs), ncol = n )
prec_weight <- matrix(0, nrow = length(probs), ncol = n )

## For every prediction day ...
for( i in 1:n )
{
    ## Extract uncertainty of all predictors and find its quantiles
    v <- var_table[ , i]
    Q <- quantile( v, probs )

    ## Counter. Counts which quantile are we in.
    j <- 1

    ## For every quantile ...
    for( q in Q  )
    {
        ## Find who is in that quantile, get corresponding
        ## predictions and get corresponding weights, both
        ## exponential and precision.
        ind <- which( v < q )
        p <- pred_table[ ind, i]
        exp_w  <- exp(-v[ind])
        prec_w <- 1/v[ind] 

        ## Weight the predictors according to the weights and store
        ## the weighted prediction in the preallocated matrices.
        exp_weight [j,i] <- sum(p* exp_w) / sum( exp_w)
        prec_weight[j,i] <- sum(p*prec_w) / sum(prec_w)

        ## Increment counter and repeat.
        j <- j + 1
    }  
}

## Prealllocate memory - every quantile has a corresponding prediction
## skill
exp_rhos  <- numeric(length(probs))
prec_rhos <- numeric(length(probs))

## For every quantile...
for( j in 1:length(probs) )
{
    ## Find the skill of the predictor built from that quantile in
    ## previous loop.
    exp_rho  <- cor( exp_weight[j, ], chl_p1wk, use = "complete.obs")
    prec_rho <- cor(prec_weight[j, ], chl_p1wk, use = "complete.obs")

    ## Store in memory.
    exp_rhos[j]  <- exp_rho
    prec_rhos[j] <- prec_rho
}


pdf(paste0("plots/exp_weighted_predictions_", tolower(xtension), ".pdf") )
plot(probs,
     exp_rhos,
     type = "l",
     xlab = "Quantile",
     ylab = "Rho",
     main = paste0( xtension, " skill per quantile for exp-weights" )
     )
dev.off()

pdf(paste0("plots/precision_weighted_predictions_", tolower(xtension), ".pdf") )
plot(probs,
     prec_rhos,
     type = "l",
     xlab = "Quantile",
     ylab = "Rho",
     main = paste0( xtension, " skill per quantile for precision-weights" )
     )
dev.off()

## First row is exp weight, second row is precision weight
predictions <- matrix( 0, nrow = 2, ncol = n )
cutoffs <- c(
    probs[which.max(exp_rhos)],
    probs[which.max(prec_rhos)]
)

## j == 1 is exp, j == 2 is precision 
for( j in 1:2 )
{
    
    ## j == 1 is exp weights, j == 2 is precision weights
    if( j == 1 )
        weight_func  <- function(x) exp(-x)
    else
        weight_func <- function(x) (1 / x)
        
    ## For every prediction day ...
    for( i in 1:n )
    {

        ## Extract uncertainty of all predictors and find its quantiles
        v <- var_table[ , i]
        q  <- quantile( v, cutoffs[j] )

        ## Find who is in that quantile, get corresponding
        ## predictions and get corresponding weights, both
        ## exponential and precision.
        ind <- which( v < q )
        p <- pred_table[ ind, i]

        w  <- weight_func(v[ind])
        
        ## Weight the predictors according to the weights and store
        ## the weighted prediction in the preallocated matrices.
        predictions[j,i] <- sum(p * w) / sum(w)
    }
}

pdf( paste0("plots/exp_time_series_", tolower(xtension), ".pdf") )
plot(chl_p1wk,
     type='l', 
     bty = "l",
     col = "black",
     ylim = c(-1, 10)
     )
lines(predictions[1, ],
      col = "red")
abline(norm_threshold,
       0,
       col = "green" )
legend("topleft",
       legend = c( "True Chl", "Exp weighted", "95% threshold" ),
       lty = 1,
       bty = 'n',
       cex = 1.25,
       col = c( "black", "red", "green" ) )
dev.off()

pdf( paste0("plots/precision_time_series_", tolower(xtension), ".pdf") )
plot(chl_p1wk,
     type='l', 
     bty = "l",
     col = "black",
     ylim = c(-1, 10)
     )
lines(predictions[2, ],
      col = "red")
abline(norm_threshold,
       0,
       col = "green" )
legend("left",
       legend = c( "True Chl", "Precision weighted", "95% threshold" ),
       lty = 1,
       bty = 'n',
       cex = 1.25,
       col = c( "black", "red", "green" ) )
dev.off()
