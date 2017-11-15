#!/usr/bin/Rscript
library(rEDM)
library(zoo)
source( "../helpers/helper.r" )

## Get a weighted prediction predictors and their uncertainty
mve_pred <- function(pred_table,
                rhos)
{
    ord <- order( rhos, decreasing = TRUE )
    ind <- ord[ 1 : ceiling(sqrt(length(rhos))) ]
    ret <- colMeans( pred_table[ ind, ] )
    return( ret )
}

## Get a uncertainty-weighted prediction
uwe_pred <- function(pred_table, var_table)
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


## If you want to download the file for some reason... don't!
## filename = paste0("http://science.sciencemag.org/highwire/filestream/683325/",
##                   "field_highwire_adjunct_files/1/aag0863_SupportingFile_Suppl._Excel_seq1_v2.xlsx")

## Load data and immediately normalize to zero mean and unit std dev
df <-normalize_df(read.csv("originals/data.csv",
                           header = TRUE,
                           sep = "," )
                  )

## Parameters
n_vars <- ncol(df) ## == 3 cuz x, y, z
E <- 3
n_lags <- 3 ## Lags: 0, -1, -2
n_samp <- 3 ##100 ## Number of random libraries
pred <- c(2501,3000) ## Prediciton always constant
pred_size <- pred[2] - pred[1] + 1
lib_sizes <- (1:3)*10

## Library changes in size and is also random
uwe_results <- data.frame(lib_sizes = lib_sizes)
mve_results <- data.frame(lib_sizes = lib_sizes)

## Craeate lags for every variable and make y_p1
df <- lag_every_variable(df, n_lags)
y_p1 <- c( as.vector(lag(zoo(df$y),1)) , NA )
df <- data.frame(y_p1 = y_p1, df )


## We know the allowed combinations must have one unlagged
## variable. By the way function combn is structured, these will be
## the first n_comb (defined below).
combinations <- combn( n_vars*n_lags, E ) ## Get ALL cobminations.
## Only this many are OK.
n_comb <- choose( n_vars*n_lags, E )- choose(n_vars*(n_lags-1), E ) 
## Throw away the rest.
combinations <- combinations[ ,1:n_comb ]

## Preallocate memory for predictions and uncertainties
pred_table <- matrix( NA,  nrow = n_comb, ncol = pred_size )
var_table  <- matrix( Inf, nrow = n_comb, ncol = pred_size )
rhos       <- matrix( NA,  nrow = 2,      ncol = n_samp    )
prediction <- numeric(                           pred_size )
ranks      <- numeric(            n_comb                   )

for (lib_ind in 1:nrow(results) )
{
    ## For every random library starting point
    for( smp in 1:n_samp )
    {
        ## Choose random lib of length exactly 100. Pred is always the same
        ## and they never overlap.
        lib  <- sample(501:2001, 1)
        lib  <- c(lib, lib + results$lib_size[lib_ind]- 1 )
        
        for( i in 1:n_comb )
        {
            cols <- combinations[,i]
            
            ## Just to make sure, all cols have one unlagged coordinate.
            stopifnot( min(cols) <= n_vars )
            cols <- cols + 1 ## cuz first col is y_p1
            
            output <- block_lnlp(df,
                                 lib = lib,   ## lib is chosen randomly above
                                 pred = pred, ## pred is ALWAYS the same
                                 method = "simplex",
                                 tp = 0, ## Zero step cuz it is built into y_p1
                                 columns = cols, 
                                 target_column = 1, ## == y_p1
                                 first_column_time = FALSE,
                                 stats_only = FALSE )
            
            pred_table[i, ] <- output[[1]]$model_output$pred
            var_table [i, ] <- output[[1]]$model_output$pred_var

            ranks[i] <- block_lnlp(df,
                                   lib = lib,   ## lib is chosen randomly above
                                   pred = lib, ## pred is ALWAYS the same
                                   method = "simplex",
                                   tp = 0, ## Zero step cuz it is built into y_p1
                                   columns = cols, 
                                   target_column = 1, ## == y_p1
                                   first_column_time = FALSE,
                                   ##short_output = TRUE,
                                   stats_only = TRUE )$rho
       
        }
        
        truth <- output[[1]]$model_output$obs
        uwe <- uwe_pred( pred_table, var_table)
        mve <- mve_pred( pred_table, ranks ) 
        rhos[1, smp] <- cor(uwe, truth, use = "complete.obs" )
        rhos[2, smp] <- cor(mve, truth, use = "complete.obs" )
    }

    quarts <- quantile(rhos[1, ], probs = c(0.25,0.5,0.75))
    uwe_results$avg[lib_ind] <- mean(rhos[1, ])
    uwe_results$bot[lib_ind] <- quarts[1]
    uwe_results$med[lib_ind] <- quarts[2]
    uwe_results$top[lib_ind] <- quarts[3]

    quarts <- quantile(rhos[2 ], probs = c(0.25,0.5,0.75))
    mve_results$avg[lib_ind] <- mean(rhos[2, ])
    mve_results$bot[lib_ind] <- quarts[1]
    mve_results$med[lib_ind] <- quarts[2]
    mve_results$top[lib_ind] <- quarts[3]

    
}

write.csv(uwe_results,
          file = "data/uwe_rhos.csv")

write.csv(mve_results,
          file = "data/mve_rhos.csv")


## noo <- read.csv("data/rhos.csv", header = TRUE, sep = "," )
## print(noo)
