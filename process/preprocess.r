#!/usr/bin/Rscript
cent <- function( x )
    return( mean( x, na.rm = TRUE ) )
scal <- function( x )
    return( sd(x, na.rm = TRUE ) )
    
## Load Melissa's raw data from file
## Prepare the data.
load( "originals/chl_block_full.Rdata" )

## Create a data frame that will contain the causal variables we want
## to consider.
df           <- data.frame( serial_day = orig_block$serial_day )

chl_c <- cent(   orig_block$chlA )
chl_s <- scal(   orig_block$chlA )
df$chl_p1wk <- ( orig_block$chlA_.1wk - chl_c ) / chl_s
df$chl      <- ( orig_block$chlA      - chl_c ) / chl_s
df$chl_1wk  <- ( orig_block$chlA_1wk  - chl_c ) / chl_s
df$chl_2wk  <- ( orig_block$chlA_2wk  - chl_c ) / chl_s

silicate_c <- cent(  orig_block$silicate )
silicate_s <- scal(  orig_block$silicate )
df$silicate     <- ( orig_block$silicate     - silicate_c ) / silicate_s
df$silicate_1wk <- ( orig_block$silicate_1wk - silicate_c ) / silicate_s
df$silicate_2wk <- ( orig_block$silicate_2wk - silicate_c ) / silicate_s

nitrate_c <- cent(  orig_block$nitrate )
nitrate_s <- scal(  orig_block$nitrate )
df$nitrate     <- ( orig_block$nitrate      - nitrate_c ) / nitrate_s
df$nitrate_1wk <- ( orig_block$nitrate_1wk  - nitrate_c ) / nitrate_s
df$nitrate_2wk <- ( orig_block$nitrate_2wk  - nitrate_c ) / nitrate_s

nitrite_c <- cent(  orig_block$nitrite )
nitrite_s <- scal(  orig_block$nitrite )
df$nitrite     <- ( orig_block$nitrite      - nitrite_c ) / nitrite_s
df$nitrite_1wk <- ( orig_block$nitrite_1wk  - nitrite_c ) / nitrite_s
df$nitrite_2wk <- ( orig_block$nitrite_2wk  - nitrite_c ) / nitrite_s

AvgTemp_1wk_c <- cent(  orig_block$AvgTemp_1wk )
AvgTemp_1wk_s <- scal(  orig_block$AvgTemp_1wk )
df$AvgTemp_1wk     <- ( orig_block$AvgTemp_1wk     - AvgTemp_1wk_c ) / AvgTemp_1wk_s
df$AvgTemp_1wk_1wk <- ( orig_block$AvgTemp_1wk_1wk - AvgTemp_1wk_c ) / AvgTemp_1wk_s
df$AvgTemp_1wk_2wk <- ( orig_block$AvgTemp_1wk_2wk - AvgTemp_1wk_c ) / AvgTemp_1wk_s

AvgDens_1wk_c <- cent( orig_block$AvgDens_1wk )
AvgDens_1wk_s <- scal( orig_block$AvgDens_1wk )
df$AvgDens_1wk     <- (orig_block$AvgDens_1wk     - AvgDens_1wk_c ) / AvgDens_1wk_s
df$AvgDens_1wk_1wk <- (orig_block$AvgDens_1wk_1wk - AvgDens_1wk_c ) / AvgDens_1wk_s
df$AvgDens_1wk_2wk <- (orig_block$AvgDens_1wk_2wk - AvgDens_1wk_c ) / AvgDens_1wk_s

U_WIND_c <- cent( orig_block$U_WIND )
U_WIND_s <- scal( orig_block$U_WIND )
df$U_WIND     <- (orig_block$U_WIND     - U_WIND_c ) / U_WIND_s
df$U_WIND_1wk <- (orig_block$U_WIND_1wk - U_WIND_c ) / U_WIND_s
df$U_WIND_2wk <- (orig_block$U_WIND_2wk - U_WIND_c ) / U_WIND_s

## Reorder the data frame
df <- df[ order(df$serial_day), ]
row.names( df ) <- 1:nrow( df )

## Make sure serial_day increases in time
stopifnot( all(
    diff( df$serial_day ) > 0
) )

save(df,
     norm_threshold,
     file = paste0("data/processed_block.Rdata" )
     )

