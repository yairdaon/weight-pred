#!/usr/bin/Rscript
library(rEDM)
source( "helpers/helper.r" )
source( "helpers/plotting.r" )

## Load the data: df, chl_threshold and norm_threshold.
load( "data/processed_block.Rdata" )

args <- commandArgs(trailingOnly = TRUE)
if( length(args) == 0 || grepl("bloom", args, ignore.case = TRUE )  ) {
    df <- get_blooms( df, norm_threshold, c(min(df$serial_day), max(df$serial_day) ) )
    xtension <- "bloom"
} else {
    xtension <- "full"
}

main_output <- block_lnlp(df,
                          method = "s-map",
                          tp = 0,             
                          columns = c("chl", "silicate_1wk", "AvgDens", "silicate" ),
                          target_column = "chl_p1wk",
                          theta =  8, 
                          num_neighbors = -1,
                          stats_only = FALSE )

serial_day <- df$serial_day
main_pred <- main_output[[1]]$model_output$pred
stopifnot(
    length(main_pred) == length(serial_day)
)

## The interpolator 
pred_fun <- approxfun(serial_day, main_pred) 

## The data, after subtractiing the chlorophyll data
leftover          <- df
leftover$chl_p1wk <- leftover$chl_p1wk - pred_fun( leftover$serial_day      )
leftover$chl      <- leftover$chl      - pred_fun( leftover$serial_day - 7  )
leftover$chl_1wk  <- leftover$chl_1wk  - pred_fun( leftover$serial_day - 14 )
leftover$chl_2wk  <- leftover$chl_2wk  - pred_fun( leftover$serial_day - 21 )

mu  <- mean(leftover$chl_p1wk, na.rm = TRUE)
sig <- sd  (leftover$chl_p1wk, na.rm = TRUE) 

leftover$chl_p1wk <- (leftover$chl_p1wk - mu) / sig
leftover$chl      <- (leftover$chl      - mu) / sig
leftover$chl_1wk  <- (leftover$chl_1wk  - mu) / sig
leftover$chl_2wk  <- (leftover$chl_2wk  - mu) / sig

## Start at 4, so we skip serial_day, chl_p1wk and chl
other_vars <- names(df)[4:length(names(df))]
n_vars <- length( other_vars )
combinations <- combn( n_vars, 3 )
leftover_pred <- matrix(NA, nrow=ncol(combinations), ncol=nrow(df) )

for( i in 1:ncol(combinations) )
{
    comb <- combinations[,i]
    ## cols <- c( "chl", combinations[,i] )
    cols <- c( 3, comb + 3 )
      
    output <- block_lnlp(leftover,
                         method = "s-map",
                         tp = 0, # Time shift is built into the target
                         columns = cols, ##c( 3, comb + 3 ),
                         target_column = 2, ##"chl_p1wk",
                         theta = 8,
                         num_neighbors = -1, ## Take all neighbours
                         first_column_time = FALSE,
                         short_output = FALSE,
                         stats_only = FALSE )

    ## Rescale the predictions so they are with the same unit as
    ## previous data.
    leftover_pred[i, ] <- output[[1]]$model_output$pred * sig + mu
}

save(leftover_pred,
     main_pred,
     other_vars,
     combinations,
     file = paste0("data/leftover_predictions_", xtension, "_loocv.Rdata" ) )
