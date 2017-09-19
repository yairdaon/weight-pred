#!/usr/bin/env Rscript
df <- data.frame()
for(year in 1991:1993 )
{

    file_df <- read.table(paste0("http://rredc.nrel.gov/solar/old_data/nsrdb/1991-2010/data/hourly/722900/722900_", year, "_solar.csv"),
                          sep = ",",
                          header = TRUE )
    df <- rbind( df, file_df )

}
times <- paste0( df$YYYY.MM.DD, " ", df$HH.MM..LST )
times <- as.POSIXct( times )
