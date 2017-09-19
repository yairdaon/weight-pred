#!/usr/bin/Rscript
cent <- function( x )
    return( mean( x, na.rm = TRUE ) )
scal <- function( x )
    return( sd(x, na.rm = TRUE ) )
normalize <- function( x, mu, sig )
    return( (x-mu)/sig )

## Load Hao's data from file. Add my solar radiation data to the same
## data frame. 
load( "originals/chl_block_full.Rdata" )
load( "data/solar_radiation.Rdata" )
orig_block <- merge(x = orig_block, y = rad_block, by = "serial_day", all.x = TRUE)

## Create a data frame that will contain the causal variables we want
## to consider.
df <- data.frame( serial_day = orig_block$serial_day )

c <- cent(   orig_block$chlA )
s <- scal(   orig_block$chlA )
df$chl_p1wk <- ( orig_block$chlA_.1wk - c ) / s
df$chl      <- ( orig_block$chlA      - c ) / s
df$chl_1wk  <- ( orig_block$chlA_1wk  - c ) / s
df$chl_2wk  <- ( orig_block$chlA_2wk  - c ) / s

c <- cent(  orig_block$silicate )
s <- scal(  orig_block$silicate )
df$silicate     <- ( orig_block$silicate     - c ) / s
df$silicate_1wk <- ( orig_block$silicate_1wk - c ) / s
df$silicate_2wk <- ( orig_block$silicate_2wk - c ) / s

c <- cent(  orig_block$nitrate )
s <- scal(  orig_block$nitrate )
df$nitrate     <- ( orig_block$nitrate      - c ) / s
df$nitrate_1wk <- ( orig_block$nitrate_1wk  - c ) / s
df$nitrate_2wk <- ( orig_block$nitrate_2wk  - c ) / s

c <- cent(  orig_block$nitrite )
s <- scal(  orig_block$nitrite )
df$nitrite     <- ( orig_block$nitrite      - c ) / s
df$nitrite_1wk <- ( orig_block$nitrite_1wk  - c ) / s
df$nitrite_2wk <- ( orig_block$nitrite_2wk  - c ) / s

c <- cent(  orig_block$AvgTemp_1wk )
s <- scal(  orig_block$AvgTemp_1wk )
df$AvgTemp     <- ( orig_block$AvgTemp_1wk     - c ) / s
df$AvgTemp_1wk <- ( orig_block$AvgTemp_1wk_1wk - c ) / s
df$AvgTemp_2wk <- ( orig_block$AvgTemp_1wk_2wk - c ) / s

c <- cent( orig_block$AvgDens_1wk )
s <- scal( orig_block$AvgDens_1wk )
df$AvgDens     <- (orig_block$AvgDens_1wk     - c ) / s
df$AvgDens_1wk <- (orig_block$AvgDens_1wk_1wk - c ) / s
df$AvgDens_2wk <- (orig_block$AvgDens_1wk_2wk - c ) / s

c <- cent( orig_block$U_WIND )
s <- scal( orig_block$U_WIND )
df$U_WIND     <- (orig_block$U_WIND     - c ) / s
df$U_WIND_1wk <- (orig_block$U_WIND_1wk - c ) / s
df$U_WIND_2wk <- (orig_block$U_WIND_2wk - c ) / s

## c <- cent( orig_block$SolRad )
## s <- scal( orig_block$SolRad )
## df$SolRad     <- (orig_block$SolRad     - c ) / s
## df$SolRad_1wk <- (orig_block$SolRad_1wk - c ) / s
## df$SolRad_2wk <- (orig_block$SolRad_2wk - c ) / s

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

