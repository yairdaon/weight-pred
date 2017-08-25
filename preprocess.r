#!/usr/bin/Rscript
source( "helper.r" )


## Load Melissa's raw data from file
raw <- read.table("originals/data_20111121.txt",
                 header = TRUE,
                 sep = "\t",
                 ## nrows = 100,
                 na.strings ="NaN" )

## ## Print the variables' names, if you'd like to.
## names( raw )

## Fix the date and serial day issues according to some (potentially
## wrong) recipe found online. This, however, seems to agree with
## Hao's data.
dates      <- as.POSIXlt(as.Date(paste(raw$Year, raw$Month, raw$Day, sep="-"))) 
serial_day <- as.numeric(dates)/86400 + 719529

## Create a raw data frame that will contain the variables I want to
## consider.  
df           <- data.frame( dates, serial_day )

## Organisms
df$chl       <- transform( raw$Chlorophyll..surface )
df$chl_p2    <- chl_week_modifier( df$chl )
df$chl_p1    <- chl_half_modifier( df$chl )
df$dino      <- transform( raw$All.Dinoflagellates..cell.counts ) 
df$diatoms   <- transform( raw$All.Diatoms..cell.counts )

## Environmental variables:

## Incoming Water
## df$rain      <- transform( past_week(raw$Lindberg.Field.Rain) )
df$river     <- transform( past_week(raw$Los.Penasquitos.River.Flow) )

## Nutrients
df$nitrate   <- transform( past_week(raw$Nitrate..surface,   raw$Nitrate..bottom   ) )
df$phosphate <- transform( past_week(raw$Phosphate..surface, raw$Phosphate..bottom ) )
df$silicate  <- transform( past_week(raw$Silicate..surface,  raw$Silicate..surface ) )
df$nitrite   <- transform( past_week(raw$Nitrite..surface,   raw$Nitrite..surface  ) )
df$ammonia   <- transform( past_week(raw$Ammonia..surface,   raw$Ammonia..surface  ) )   
## df$nitro     <- transform( df$nitrate + df$nitrite + df$ammonia ) ONLY NAs!!! 
df$NO_total  <- transform( df$nitrate + df$nitrite )
df$NO_spread <- transform( df$nitrate - df$nitrite )

## Wind
df$wind_v    <- transform( past_week(raw$V.wind.componet.at.10m) )
df$wind_u    <- transform( past_week(raw$U.wind.componet.at.10m) )

## Sea properties
df$density   <- transform( past_week(raw$Surface.density.using.Shore.Station.Program.T..S,             
                                       raw$Bottom.density.using.Shore.Station.Program.T.S) )
df$temp      <- transform( past_week(raw$Daily.surface.water.temperature,
                                       raw$Daily.bottom.water.temperature) )
df$salinity  <- transform( past_week(raw$Daily.surface.salinity,
                                       raw$Daily.bottom.salinity) )

write.csv(df, file = "data/raw.csv",
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

## A boolean vector with TRUE on days we have chlorophyll data
indices <- !is.na(df$chl)

## A vector of MATLAB serial days where we have chlorophyll data
chl_days <- df[ !is.na(df$chl), ]$serial_day

## Day differences between consecutive chlorophyll observations.
step <- c( 3, diff( chl_days ) )

## The indices where the time difference is NOT 3 or our days
indices <- unique( c( which( step < 3 ), which( step > 4 ) ) )

## days that should be an NA row, since the day after, there is a
## measurement that is not within 3,4 days of previous measurement
NA_row_days <- chl_days[indices] - 1

## Restrict the dataframe to the rows we need
df <- df[ df$serial_day %in% c( chl_days, NA_row_days ), ]
df$step <- c( 3, diff( df$serial_day ) )

## Make sure we insert those NA rows!!!
df[ df$serial_day %in% NA_row_days, ] <- rep( NA, ncol(df) )


write.csv( df, file = "data/univariate_data.csv",
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


