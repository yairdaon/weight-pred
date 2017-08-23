#!/usr/bin/Rscript
library(rEDM)
source( "helper.r" )
    
predict <- function(past, ## What we use for prediction 
                    future, ## The thing we try to predict
                    E,
                    tp,
                    which_obs = "all" )
{
    ## Lag it... 
    if( past == "chl" )
        past <- multi_lag_chl( E )
    else
        past <- multi_lag_var(past,
                              E )
    
    ## ... get the column names for later use ...
    pred_cols <- names( past )
    pred_cols <- pred_cols[-1]

    ## ## ... and remove all illegal indices.
                 ## indices <- legal_indices(restricted_df$serial_day,
    ##                          E,
    ##                          which_obs )
    ## lagged_predictor <- lagged_predictor[ indices, ]
    
    ## Predictee:
    
    ## Get its forwarded data frame

    if( future == "chl" && tp %% 2 == 0 )
        future <- half_week_steps( week_chl, tp )
    
    
    
    df <- merge(x = past, y = future, by = "serial_day", all.x = TRUE)
    df <- na.omit( df ) 
    ## df <- restrict( df )
    print( df[1:12, ] )
    output <- block_lnlp(df,
                         method  = "simplex",
                         columns = pred_cols,
                         tp = 0,
                         target_column = 1,
                         stats_only = TRUE,
                         first_column_time = TRUE,
                         silent = FALSE)

    
    return( output$rho )
           
}

predictor <- "chl"
predictee <- "chl"
E         <- 4
tp        <- 0
which_obs <- "all"

rhos <- list()
for( e in E )
{
    rho <- predict(predictor,
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
