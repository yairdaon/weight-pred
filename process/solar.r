#!/usr/bin/env Rscript
library( zoo )
source( "helpers/helper.r" )

df <- data.frame()
for( year in 1991:2010 )
{
    new_df <- read.table(paste0("http://rredc.nrel.gov/solar/old_data/nsrdb/1991-2010/data/hourly/722900/722900_", year, "_solar.csv"),
                     sep = ",",
                     header = TRUE )
    df <- rbind( df, new_df )
}

## Create the MATLAB serial day time stamp
times <- paste0( df$YYYY.MM.DD, " ", df$HH.MM..LST )
times <- as.POSIXct( times )
times <- round( times, "day" )
df$serial_day <-date2serial( times )

## Sum solar radiation by day
df <- aggregate(df$METSTAT.Glo..Wh.m.2.,
                       by=list(serial_day=df$serial_day),
                       FUN=sum )
names(df)[names(df) == 'x'] <- 'daily'
stopifnot( all(
    df$daily > 0
) )

## Get sum of total radiation over past week
df$weekly <- rollapply(df$daily,
                       width = 7,
                       FUN = sum,
                       na.rm = TRUE,
                       fill = NA,
                       align = "right" )

rad_block <- data.frame( serial_day = df$serial_day )
ts <- zoo( df$weekly )
rad_block$SolRad     <-      ts
rad_block$SolRad_1wk <- lag( ts,  -7,  na.pad = TRUE )
rad_block$SolRad_2wk <- lag( ts, -14, na.pad = TRUE )

stopifnot(all(
        diff( rad_block$serial_day ) == 1
) )

save(rad_block,
     file = "data/solar_radiation.Rdata" )

