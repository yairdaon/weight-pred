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

all_pred_obs  <- df$chl_p1wk[pred[1]:pred[2]]
all_pred_time <- pred[1]:pred[2]
T             <- all_pred_time[1] - 1
n             <- length( all_pred_time )

vars <- c("nitrate",
          "nitrate_m1wk",
          "nitrate_m2wk",
          "silicate",
          "silicate_m1wk",
          "silicate_m2wk",
          "nitrite",
          "nitrite_m1wk",
          "nitrite_m2wk",
          "AvgTemp_1wk",
          "AvgTemp_1wk_m1wk",
          "AvgTemp_1wk_m2wk",
          "AvgDens_1wk",
          "AvgDens_1wk_m1wk",
          "AvgDens_1wk_m2wk",
          "U_WIND",
          "U_WIND_m1wk",
          "U_WIND_m2wk",
          "chl",
          "chl_m1wk",
          "chl_m2wk" )


best <- block_lnlp(df,
                   lib = lib,
                   pred = pred,
                   method = "simplex", ##"s-map",
                   tp = 0, # Cuz the time shift is alread in the
                                        # target
                   columns = c("chl",
                               "silicate_m1wk",
                               "AvgDens_1wk_m1wk",
                               "silicate"),
                   target_column = "chl_p1wk",
                   #theta = 1.8,
                   short_output = TRUE,
                   stats_only = FALSE )


best_pred       <- numeric(n)
time            <- best[[1]]$model_output$time - T
best_pred[time] <- best[[1]]$model_output$pred


min_vars   <- rep( Inf, n )
var_pred   <- rep( NA,  n )
cum_pred   <- numeric(n)
cum_weight <- numeric(n)

## Should have |vars|/4!4! = 8!/4!4! = ... = 70 combinations
combinations <- combn( vars, 4 )

for( i in 1:dim(combinations)[2] )
{
    cols <- combinations[,i]

    output <- block_lnlp(df,
                         lib = lib,
                         pred = pred,
                         method = "simplex", ##"s-map",
                         tp = 0, # Cuz the time shift is alread in the
                                 # target
                         columns = cols,
                         target_column = "chl_p1wk",
                         ##theta = 1.5,
                         short_output = TRUE,
                         stats_only = FALSE )

    time <- output[[1]]$model_output$time - T

    ## The model's predictions
    pr <- numeric(n)
    pr[time] <- output[[1]]$model_output$pred

    ## The prediction's associated variances
    vars <- rep( Inf, n )
    vars[time] <- output[[1]]$model_output$pred_var

    ## Sanity check
    stopifnot( all(
    (vars==Inf) == (pr==0)
    ) )

    ## Choose new predictions if they are less uncertain than previous ones
    min_ind           <- which( vars < min_vars )
    min_vars[min_ind] <- vars[min_ind]
    var_pred[min_ind] <- pr[min_ind]

    ## Weigh new predictions compare to previous ones
    weight     <- ceiling(exp(-vars/1000))
    cum_pred   <- cum_pred + weight * pr 
    cum_weight <- cum_weight + weight
    
}

cum_pred <- cum_pred / cum_weight


## x11()
## plot(all_pred_time,
##      all_pred_obs,
##      type = "l",
##      col = "black")
## lines(prediction, col = "red" )
## lines(exp_prediction, col = "green" )
## lines(best_pred, col = "blue" )
## hold(1)

print( paste0( "Out of sample rho using best 4D model: "         ,         cor(best_pred, all_pred_obs, use = "complete.obs") ) )
print( paste0( "Out of sample rho using weighted 4D predictors: ",         cor(cum_pred,  all_pred_obs, use = "complete.obs") ) )
print( paste0( "Out of sample rho using minimal variance 4D predictors: ", cor(var_pred,  all_pred_obs, use = "complete.obs") ) )

warnings()
