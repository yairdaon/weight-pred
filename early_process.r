#!/usr/bin/Rscript
source( "helper.r" )


## Load Melissa's raw data from file
raw <- read.table("originals/data_20111121.txt",
                 header = TRUE,
                 sep = "\t",
                 ## nrows = 100,
                 na.strings ="NaN" )

## Print the variables names, if you want
## names( raw )

## Fix the date and serial day issues according to some (potentially
## wrong) recipe found online. This, however, seems to agree with
## Hao's data.
dates      <- as.POSIXlt(as.Date(paste(raw$Year, raw$Month, raw$Day, sep="-"))) 
serial_day <- as.numeric(dates)/86400 + 719529

## Create a raw data frame that will contain the variables I want to
## consider.  
mine           <- data.frame( dates, serial_day )

## Organisms
mine$chl       <- transform( raw$Chlorophyll..surface )
mine$chl_week  <- chl_week_modifier( mine$chl )
mine$chl_half  <- chl_half_modifier( mine$chl )
mine$dino      <- transform( raw$All.Dinoflagellates..cell.counts ) 
mine$diatoms   <- transform( raw$All.Diatoms..cell.counts )

## Environmental variables:

## Incoming Water
## mine$rain      <- transform( past_week(raw$Lindberg.Field.Rain) )
mine$river     <- transform( past_week(raw$Los.Penasquitos.River.Flow) )

## Nutrients
mine$nitrate   <- transform( past_week(raw$Nitrate..surface,   raw$Nitrate..bottom   ) )
mine$phosphate <- transform( past_week(raw$Phosphate..surface, raw$Phosphate..bottom ) )
mine$silicate  <- transform( past_week(raw$Silicate..surface,  raw$Silicate..surface ) )
mine$nitrite   <- transform( past_week(raw$Nitrite..surface,   raw$Nitrite..surface  ) )
mine$ammonia   <- transform( past_week(raw$Ammonia..surface,   raw$Ammonia..surface  ) )   
## mine$nitro     <- transform( mine$nitrate + mine$nitrite + mine$ammonia ) ONLY NAs!!! 
mine$NO_total  <- transform( mine$nitrate + mine$nitrite )
mine$NO_spread <- transform( mine$nitrate - mine$nitrite )

## Wind
mine$wind_v    <- transform( past_week(raw$V.wind.componet.at.10m) )
mine$wind_u    <- transform( past_week(raw$U.wind.componet.at.10m) )

## Sea properties
mine$density   <- transform( past_week(raw$Surface.density.using.Shore.Station.Program.T..S,             
                                       raw$Bottom.density.using.Shore.Station.Program.T.S) )
mine$temp      <- transform( past_week(raw$Daily.surface.water.temperature,
                                       raw$Daily.bottom.water.temperature) )
mine$salinity  <- transform( past_week(raw$Daily.surface.salinity,
                                       raw$Daily.bottom.salinity) )

## Now we show our transformed data-set, before choosing Hao's rows
for( var in names( mine ) )
    if( FALSE ) ## PlotOrNot?
        harry_plotter( mine[ ,var ], var )

write.csv( mine, file = "data/processed_daily_data.csv", quote = FALSE )

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
