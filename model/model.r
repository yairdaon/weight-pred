#!/usr/bin/Rscript
library(rEDM)
library(zoo)
source( "../helpers/helper.r" )


## In this script we show prediction skill of the uncertainty -
## weighting scheme. We randomly sample libraries from synthetic data
## and present their skill at out-of-sample prediction.

## If you want to download the file for some reason... don't!
## filename = paste0("http://science.sciencemag.org/highwire/filestream/683325/",
##                   "field_highwire_adjunct_files/1/aag0863_SupportingFile_Suppl._Excel_seq1_v2.xlsx")

predict <- function(raw_df,
                    variables = names(raw_df), 
                    E = 3, ## Embedding dimension of the system.
                    n_lags = 3, ## 0, -1,..., -(n_lags-1)
                    n_samp = 5, ## Number of random libraries, should be in the hundreds
                    lib_range = c(501:2001),  ## Library set.
                    pred = c(2501,3000), ## Prediciton set.
                    lib_sizes = (1:3)*10     ## Library changes in size and is also random
                    )
{   
    ## Immediately rescale the data frame and keep track of the
    ## normalizing factors.
    mus  <- get_mus(raw_df)
    sigs <- get_sigs(raw_df)
    df   <- data.frame(scale(raw_df))
      
    ## number of variables in the dynamical system
    n_vars <- ncol(df)

    ## Preallocate 
    results <- data.frame(lib_sizes = lib_sizes)
    results$avg <- numeric(nrow(results))
    results$top <- numeric(nrow(results))
    results$med <- numeric(nrow(results))
    results$bot <- numeric(nrow(results))

    ## Create lags for every variable and make y_p1
    df <- lag_every_variable(df, n_lags)
    

    ## Create the list of combinations
    combinations <- make_combinations(n_vars, n_lags, E)
    n_comb <- ncol(combinations)
    
    ## Preallocate memory for predictions and uncertainties
    pred_size   <- pred[2] - pred[1] + 1
    pred_table  <- matrix( NA,  nrow = n_comb, ncol = pred_size )
    var_table   <- matrix( Inf, nrow = n_comb, ncol = pred_size )
    rhos        <- rep   ( NA,         n_samp                   )
    prediction  <- numeric(                           pred_size )

    ## Do the analysis for every variable. 
    for( curr_var in variables )
    {
        for (lib_ind in 1:nrow(results) )
        {
            ## For every random library starting point
            for( smp in 1:n_samp )
            {
                ## Choose random lib
                lib <- random_lib(lib_range, lib_sizes[lib_ind])
                
                for( comb in 1:n_comb )
                {

                    ## Just to make sure, all variable sets have one
                    ## unlagged coordinate.
                    ## stopifnot( min(combinations[,comb]) <= n_vars )
                    
                    ## Variables used for prediction are chosen
                    ## according to the current combinations. Since
                    ## first cols are x_p1, y_p1, z_p1, we skip them
                    ## and add n_vars.
                    cols <- combinations[, comb] + n_vars
                    target_column <- paste0( curr_var, "_p1" )

                    output <- block_lnlp(df,
                                         lib = lib,   ## lib is chosen randomly above
                                         pred = pred, ## pred is ALWAYS the same
                                         method = "simplex",
                                         tp = 0, ## Lags built into df
                                         columns = cols, 
                                         target_column = target_column,
                                         first_column_time = FALSE,
                                         short_output = FALSE,
                                         stats_only = FALSE )

                    ## time <- output[[1]]$model_output$time
                    pr   <- output[[1]]$model_output$pred[pred[1]:pred[2]]
                    vars <- output[[1]]$model_output$pred_var[pred[1]:pred[2]]

                    ## rho <- output[[1]]$stats$rho
                    ## calc_rho <- cor(pr, df[time,var_ind], use = "complete.obs" ) 
                    ## Not sure why these are not exactly equal...
                    ## if(  abs(calc_rho-rho) > 1e-14 )
                    ##     stop( calc_rho - rho )

                    ## Save ALL data.
                    pred_table[comb, ] <- pr 
                    var_table [comb, ] <- vars
                }


                mu  <- mus$curr_var
                sig <- sigs$curr_var
                
                prediction <- weighted_prediction(var_table, pred_table)
                truth <- df[ pred[1]:pred[2], curr_var ]

                ## We don't really need to shift and rescale variables
                ## back, since the correlation is invariant to shifts
                ## and scaling (when operations performed on both
                ## variables).
                prediction <- descale(prediction, mu, sig)
                truth <- descale(truth, mu, sig)

                ## Store the correlation coefficient
                rhos[smp]  <- cor(prediction, truth, use = "complete.obs" )

            }

            quarts <- quantile(rhos, probs = c(0.25,0.5,0.75))

            results$avg[lib_ind] <- mean(rhos)
            results$bot[lib_ind] <- quarts[1]
            results$med[lib_ind] <- quarts[2]
            results$top[lib_ind] <- quarts[3]
        }

        write.csv(results,
                  file = paste0("data/rhos_", curr_var, ".csv") )

        pdf( paste0("plots/rhos_", curr_var, ".pdf") )
        plot(lib_sizes,
             results$avg,
             type = "l",
             xlab = "Library Size",
             ylab = "Prediction Skill")

        lines(results$avg, col="red" )
        lines(results$avg, col="red" )
        dev.off()
        
    }
}


## Load data 
df <- read.csv("originals/data.csv",
                     header = TRUE,
                     sep = "," )


predict(df,
        variables = names(df), 
        E = 3, ## Embedding dimension of the system.
        n_lags = 3, ## 0, -1,..., -(n_lags-1)
        n_samp = 5, ## Number of random libraries, should be in the hundreds
        lib_range = c(501:2001),  ## Library set.
        pred = c(2501,3000), ## Prediciton set.
        lib_sizes = (1:10)*2 ## Library sizes
        )


## L8YTCKPF48Y070948 ## True
## L8YTCKPF98Y070945 ## Pink slip
## Outsite
