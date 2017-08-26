#!/usr/bin/Rscript
library(rEDM)
source( "plotting.r" )
source( "helper.r" )
source( "ccm.r" )

df   <- read.csv( "data/univariate_data.csv", header = TRUE )
n    <- floor(3*nrow( df )/4)

perm <- sample( nrow(df), size = nrow(df), replace = FALSE )
lib  <- perm[  1 :  n      ]
pred <- perm[(n+1):nrow(df)]

lib  <- c(   1   , n       )
pred <- c( (n+1) ,nrow(df) )

lib  <- c(1 , nrow(df) )
pred <- lib
##lib <- scan( "data/lib_days.txt" )
## lib <- which( df$serial_day %in% lib )
##print( lib )

## pred <- scan( "data/pred_days.txt" )
## pred <- which( df$serial_day %in% pred )
## print( pred )

## Univariate analysis
simplex_output <- simplex(df$chl,
                          lib = lib,
                          pred = lib,
                          tp = 1)

smap_output <- s_map(df$chl,
                     lib = lib,
                     pred = pred,
                     E = 4,
                     tp = 1)
theta <- smap_output$theta
rho <- smap_output$rho


x11()
plot(simplex_output$E,
     simplex_output$rho,
     type = "l",
     xlab = "Embedding Dimension (E)", 
     ylab = "Forecast Skill (rho)" )
hold(1)
plot(theta,
     rho,
     type = "l",
     xlab = "Nonlinearity (theta)", 
     ylab = "Forecast Skill (rho)")
hold(1)

print( paste0("Maximal rho = ", max( rho ), " with theta = ", theta[ which.max( rho ) ] ) )


warnings()
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
