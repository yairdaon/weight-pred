#!/usr/bin/Rscript
library(rEDM)
source( "plotting.r" )
source( "helper.r" )
source( "ccm.r" )

## df   <- read.csv( "originals/chl_block.csv", header = TRUE )
## n    <- floor(3*nrow( df )/4)

## perm <- sample( nrow(df), size = nrow(df), replace = FALSE )
## lib  <- perm[  1 :  n      ]
## pred <- perm[(n+1):nrow(df)]

## lib  <- c(   1   , n       )
## pred <- c( (n+1) ,nrow(df) )

## lib  <- c(1 , nrow(df) )
## pred <- lib


## ######################################
## ## Univariate analysis ###############
## ######################################
## simplex_output <- simplex(df$chl,
##                           lib = lib,
##                           pred = lib,
##                           tp = 1)

## smap_output <- s_map(df$chl,
##                      lib = lib,
##                      pred = pred,
##                      E = 4,
##                      tp = 1)
## theta <- smap_output$theta
## rho <- smap_output$rho


## x11()
## plot(simplex_output$E,
##      simplex_output$rho,
##      type = "l",
##      xlab = "Embedding Dimension (E)", 
##      ylab = "Forecast Skill (rho)" )

## hold(1)
## plot(theta,
##      rho,
##      type = "l",
##      xlab = "Nonlinearity (theta)", 
##      ylab = "Forecast Skill (rho)")
## hold(1)

## print( paste0("Maximal rho = ", max( rho ), " with theta = ", theta[ which.max( rho ) ] ) )


## warnings()


###################################
## Multivariate ###################
###################################
df   <- read.csv( "originals/chl_block.csv", header = TRUE )
lib  <- scan( "originals/lib_rows.txt"  )
pred <- scan( "originals/pred_rows.txt" )
lib  <- c( min(lib), max(lib) ) 
pred <- c( min(pred), max(pred) )

all_pred_obs  <- df$chl[pred[1]:pred[2]]
all_pred_time <- pred[1]:pred[2]

vars <- c("nitrate_m1wk",
          "nitrate_m2wk",
          "silicate_m1wk",
          "silicate",
          "nitrite_m1wk",
          "AvgTemp_1wk_m1wk",
          "AvgDens_1wk_m1wk",
          "U_WIND_m2wk",
          "chl",
          "chl_m1wk",
          "chl_m2wk" )


best <- block_lnlp(df,
                   lib = lib,
                   pred = pred,
                   method = "s-map",
                   tp = 0, # Cuz the time shift is alread in the
                                        # target
                   columns = c("chl",
                               "silicate_m1wk",
                               "AvgDens_1wk_m1wk",
                               "silicate"),
                   target_column = "chl_p1wk",
                   theta = 1.8,
                   short_output = TRUE,
                   stats_only = FALSE )


best_pred       <- numeric(length(all_pred_time))
time            <- best[[1]]$model_output$time
best_pred[time] <- best[[1]]$model_output$pred
     

cum_weight_pred <- numeric(length(all_pred_time))
cum_weight      <- numeric(length(all_pred_time))

## Should have |vars|/4!4! = 8!/4!4! = ... = 70 combinations
combinations <- combn( vars, 4 )

for( i in 1:dim(combinations)[2] )
{
    cols <- combinations[,i]
    ## print( cols )

    output <- block_lnlp(df,
                         lib = lib,
                         pred = pred,
                         method = "s-map",
                         tp = 0, # Cuz the time shift is alread in the
                                 # target
                         columns = cols,
                         target_column = "chl_p1wk",
                         theta = 1.5,
                         short_output = TRUE,
                         stats_only = FALSE )
    

    ## print( paste0("stats:", names(output[[1]]$stats)) )
    ## print( paste0("model output:", names(output[[1]]$model_output)) )
    ## print( paste0("embedding:", names(output[[1]]$embedding)) )
    ## print( paste0("params:",    names(output[[1]]$params)) )
    
    weight <- exp(-output[[1]]$model_output$pred_var)
    weight <- 1 / output[[1]]$model_output$pred_var
    weight[ is.nan(weight)  ] <- 0

    pr <- output[[1]]$model_output$pred
    pr[ is.nan(pr)  ] <- 0

    time <- output[[1]]$model_output$time
    
    cum_weight_pred[time] <- cum_weight_pred + weight * pr 

    cum_weight[time] <- cum_weight + weight
}

prediction <- cum_weight_pred / cum_weight
x11()
plot(all_pred_time,
     all_pred_obs,
     type = "l",
     col = "black")
lines(prediction, col = "red" )
lines(best_pred, col = "blue" )
hold(1)

best_rho <- 
warnings()
## lib <- scan( "data/lib_days.txt" )
## lib <- which( df$serial_day %in% lib )
## print( lib )

## pred <- scan( "data/pred_days.txt" )
## pred <- which( df$serial_day %in% pred )
## print( pred )





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
