#!/usr/bin/Rscript
date2serial <- function( date_string )
    return(
        as.numeric(as.POSIXlt(as.Date(date_string)))/86400 + 719529
    )

lib_start  <- date2serial( "2008-01-01" )
lib_end    <- date2serial( "2009-12-31" )
pred_start <- date2serial( "2010-01-01" )
pred_end   <- date2serial( "2010-06-30" )
lib_serial_days  <- c( lib_start,  lib_end)
pred_serial_days <- c(pred_start, pred_end)
save(lib_serial_days,
     pred_serial_days,     
     file="data/test_serial_days.Rdata" )

lib_start  <- date2serial( "1983-01-01" )
lib_end    <- date2serial( "2010-12-28" )
pred_start <- date2serial( "2011-01-01" )
pred_end   <- date2serial( "2012-03-30" )
lib_serial_days  <- c( lib_start,  lib_end)
pred_serial_days <- c(pred_start, pred_end)
save(lib_serial_days,
     pred_serial_days,
     file="data/oos_serial_days.Rdata" )

lib_start  <- date2serial( "1983-01-01" )
lib_end    <- date2serial( "2012-03-30" )
pred_start <- date2serial( "1983-01-01" )
pred_end   <- date2serial( "2012-03-30" )
lib_serial_days  <- c( lib_start,  lib_end)
pred_serial_days <- c(pred_start, pred_end)
save(lib_serial_days,
     pred_serial_days,
     file="data/loocv_serial_days.Rdata" )
