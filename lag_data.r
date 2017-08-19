#!/usr/bin/Rscript
library(rEDM)
source( "helper.r" )
    
predict <- function(predictor,
                    predictee,
                    E,
                    tp,
                    which_obs = "all" )
{
    file_df <- read.csv( "data/processed_data.csv", header = TRUE )

    ## Lag the predictor 
    lagged_predictor <- multi_lag_ts(file_df,
                                     predictor,
                                     E )
    
    ## Keep only legal indices
    predictor_indices <- predictor_legal_indices(df$serial_day,
                                                 E,
                                                 which_obs )
    lagged_predictor  <- lagged_predictor[ predictor_indices, ]


    predictee_ts      <- look_ahead(file_df,
                                    predictee,
                                    tp )
    
    predictee_indices <- predictee_legal_indices(df$serial_day,
                                                 tp )
    df$serial_day <- NULL
    
    pred_cols <- names( df )
        
    df[, predictee] <- file_df[[predictee]]

    ## Only keep the legal indices
    indices <- predictor_indices & predictee_indices
    
        
    output <- block_lnlp(df,
                         method  = "simplex",
                         columns = pred_cols,
                         tp = 0,
                         target_column = predictee,
                         stats_only = TRUE,
                         first_column_time = FALSE,
                         silent = TRUE)

    
    return( output$rho )
           
}

predictor <- "chl"
predictee <- "chl"
E         <- 3:7
tp        <- 0
which_obs <- "all"

rhos <- list()
for( e in E )
{
    rho <- my_ccm(predictor,
                  predictee,
                  e,
                  tp,
                  which_obs )
    
    print( paste0( predictor, " xmap ", predictee, ". Tp = ", tp, " E = ", e, " rho = ", rho ) )  
}

## x11()
## plot(E,
##      rhos,
##      type = "l",
##      col = "blue",
##      main = paste0("Prediction of ", predictee, " using ", predictor, " data using ", numLags, " lags") )
## hold(1)

## perm <- order(rho)
## print( paste0("Lags = ",
##               E[length(E)],
##               ", max rho = ",
##               rho[perm[length(rho)]],
##               ", E = ",
##               perm[length(rho)],
##               ", 2nd rho = ",
##               rho[perm[length(rho)-1]],
##               ", E = ",
##               perm[length(rho)-1] )
##       )
