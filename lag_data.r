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
    past <- multi_lag_var(past, E )
    
    ## ... get the column names for later use ...
    pred_cols <- names( past )
    pred_cols <- pred_cols[-1]
    
    ## Get its forwarded data frame
    future <- look_ahead_var( future, tp )

    target <- colnames( future )
    target <- target[-1]
    
    df <- merge(x = past, y = future, by = "serial_day", all.x = TRUE)

    ## Cleaning and restriction procedures
    ## df <- na.omit( df ) 
    df <- restrict( df )

    ## Printers
    ## print( df[1:12, ] )
    print( names( df ) )
    print( pred_cols )
    print( target )
    print( length( df$serial_day) )
    output <- block_lnlp(df,
                         method  = "simplex",
                         columns = pred_cols,
                         tp = 0,
                         target_column = target,
                         stats_only = TRUE,
                         first_column_time = TRUE,
                         silent = FALSE)

    
    return( output$rho )
           
}

predictor <- "nitrate"
predictee <- "chl"###"nitrate"
E         <- 4
tp        <- -2
which_obs <- "all"
for( E in 5 )
{
    rho <- predict(predictor,
                   predictee,
                   E,
                   tp,
                   which_obs )
    
    print( paste0( predictor, " xmap ", predictee, ". Tp = ", tp, " E = ", E, " rho = ", rho ) )  
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
