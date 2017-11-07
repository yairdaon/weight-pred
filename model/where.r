#!/usr/bin/Rscript
library(rEDM)
library(zoo)
source( "../helpers/helper.r" )
source( "../helpers/plotting.r" )

load( "runs/parameters.Rdata" )

## Get the true values we were trying to predict
raw_df <- read.csv("originals/three_species.csv",
                   header = TRUE,
                   sep = "," )
truth <- raw_df$x[pred[1]:pred[2]]
errors <- read.csv("runs/mean_tracking_errors_x.csv")

x11()
plot(pred[1]:pred[2],
     errors$weighted,
     type = "l",
     col = "blue" )
lines(pred[1]:pred[2],
      errors$mve,
      col = "red" )
hold(1)

x11()
plot(pred[1]:pred[2],
     log(errors$mve/errors$weighted),
     type = "l",
     col = "blue" )
hold(1)



## print(errors)








## libs <- read.table("runs/libraries_x.txt",
##                    header = FALSE,
##                    sep = "\t" )
