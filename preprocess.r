#!/usr/bin/Rscript
source( "helper.r" )

cent <- function( x )
    return( mean( x, na.rm = TRUE ) )
scal <- function( x )
    return( sd(x, na.rm = TRUE ) )

## Load Melissa's raw data from file
## Prepare the data.
load( "originals/chl_block_full.Rdata" )

## ## Print the variables' names, if you'd like to.
## print( names( orig_block ) )

## Create a data frame that will contain the causeal variables we want
## to consider.
df           <- data.frame( serial_day = orig_block$serial_day )

chl_c <- cent( orig_block$chlA )
chl_s <- scal(   orig_block$chlA )
df$chl_p1wk <- (orig_block$chlA_.1wk - chl_c ) / chl_s
df$chl      <- (orig_block$chlA      - chl_c ) / chl_s
df$chl_1wk  <- (orig_block$chlA_1wk  - chl_c ) / chl_s
df$chl_2wk  <- (orig_block$chlA_2wk  - chl_c ) / chl_s

silicate_c <- cent( orig_block$silicate )
silicate_s <- scal(   orig_block$silicate )
df$silicate     <- (orig_block$silicate     - silicate_c ) / silicate_s
df$silicate_1wk <- (orig_block$silicate_1wk - silicate_c ) / silicate_s
df$silicate_2wk <- (orig_block$silicate_2wk - silicate_c ) / silicate_s

nitrate_c <- cent( orig_block$nitrate )
nitrate_s <- scal(   orig_block$nitrate )
df$nitrate     <- (orig_block$nitrate      - nitrate_c ) / nitrate_s
df$nitrate_1wk <- (orig_block$nitrate_1wk  - nitrate_c ) / nitrate_s
df$nitrate_2wk <- (orig_block$nitrate_2wk  - nitrate_c ) / nitrate_s

nitrite_c <- cent( orig_block$nitrite )
nitrite_s <- scal(   orig_block$nitrite )
df$nitrite     <- (orig_block$nitrite      - nitrite_c ) / nitrite_s
df$nitrite_1wk <- (orig_block$nitrite_1wk  - nitrite_c ) / nitrite_s
df$nitrite_2wk <- (orig_block$nitrite_2wk  - nitrite_c ) / nitrite_s

AvgTemp_1wk_c <- cent( orig_block$AvgTemp_1wk )
AvgTemp_1wk_s <- scal(   orig_block$AvgTemp_1wk )
df$AvgTemp_1wk     <- (orig_block$AvgTemp_1wk     - AvgTemp_1wk_c ) / AvgTemp_1wk_s
df$AvgTemp_1wk_1wk <- (orig_block$AvgTemp_1wk_1wk - AvgTemp_1wk_c ) / AvgTemp_1wk_s
df$AvgTemp_1wk_2wk <- (orig_block$AvgTemp_1wk_2wk - AvgTemp_1wk_c ) / AvgTemp_1wk_s

AvgDens_1wk_c <- cent( orig_block$AvgDens_1wk )
AvgDens_1wk_s <- scal(   orig_block$AvgDens_1wk )
df$AvgDens_1wk     <- (orig_block$AvgDens_1wk     - AvgDens_1wk_c ) / AvgDens_1wk_s
df$AvgDens_1wk_1wk <- (orig_block$AvgDens_1wk_1wk - AvgDens_1wk_c ) / AvgDens_1wk_s
df$AvgDens_1wk_2wk <- (orig_block$AvgDens_1wk_2wk - AvgDens_1wk_c ) / AvgDens_1wk_s

U_WIND_c <- cent( orig_block$U_WIND )
U_WIND_s <- scal( orig_block$U_WIND )
df$U_WIND     <- (orig_block$U_WIND     - U_WIND_c ) / U_WIND_s
df$U_WIND_1wk <- (orig_block$U_WIND_1wk - U_WIND_c ) / U_WIND_s
df$U_WIND_2wk <- (orig_block$U_WIND_2wk - U_WIND_c ) / U_WIND_s

write.csv(df, file = "data/block.csv",
          quote = FALSE,
          row.names = FALSE )

##########################################################################
##########################################################################
##########################################################################
##########################################################################


## Hao's data.
hao <- read.csv("originals/chl_block.csv",
                header = TRUE )

## Save Hao's days to a file
write( hao$serial_day, file = "data/hao_days.txt" )

## Ceate and save the dates that correspond to Hao's
## prediction set, according to Hao's indices.
pred_rows    <- scan( "originals/pred_rows.txt" )
hao_pred     <- hao[pred_rows, ]
pred_days    <- hao_pred$serial_day
write( pred_days, file = "data/pred_days.txt" )

## Same for library
lib_rows    <- scan( "originals/lib_rows.txt" )
hao_lib     <- hao[lib_rows, ]
lib_days    <- hao_lib$serial_day
write( lib_days, file = "data/lib_days.txt" )


##########################################################################
##########################################################################
##########################################################################
##########################################################################

## Here we do the row manipulation for UNIVARIATE analysis. We make
## sure we have blocks of consecutive 3-4 day increase between
## measurements and when we don't have these, we stick a NA row, to
## denote that the block ends.

## Throw away all days where chlorophyll was not measured.
df  <- df[ !is.na(df$chl), ]

## Day differences between consecutive chlorophyll observations.
step <- c( 3, diff( df$serial_day ) )

## The indices where the time difference is NOT 3 or 4 days, sorted.
## These are the breakpoints where we add a NA row.
blocks <- sort(unique( c( which( step < 3 ), which( step > 4 ) ) ) )

## So we include data to the end of the original data frame. Just a
## programming shortcut.
blocks <- c( blocks, nrow(df)+1 )

## The NA block we separate with.
filler <- rep( NA, nrow(df) )
## filler <- matrix(data=NA,
##                  nrow=nrow(df),
##                  ncol=4 )


## The augmented data frame, with NA between observations which are
## not 3 or 4 days apart.
aug <- df[1:(blocks[1]-1), ]

for( r in 1:(length(blocks)-1) )

    ## Augment the data frame with a filler row, then the next block.
    aug <- rbind( aug, filler, df[ blocks[r]:(blocks[r+1]-1) , ] )


## Tests. Work when we fill with only one NA row.!!!

## Make sure we don't throw away ANY data
stopifnot( all(
    df$serial_day == na.omit(aug$serial_day) 
) )

## Make sure the jump before and after a NA line is not 3 or 4.
rows  <- which(is.na(aug$serial_day))
diffs <- aug$serial_day[rows+1] - aug$serial_day[rows-1]
stopifnot( all(
    (diffs > 4) | (diffs < 3) 
) )

write.csv(aug,
          file = "data/univariate_data.csv",
          quote = FALSE,
          row.names = FALSE )


##########################################################################
##########################################################################
##########################################################################
##########################################################################



## Now we show our transformed data-set, before choosing Hao's rows
for( var in names( df ) )
    if( FALSE ) ## PlotOrNot?
        harry_plotter( df[ ,var ], var )


