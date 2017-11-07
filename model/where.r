#!/usr/bin/Rscript
library(rEDM)
library(zoo)
source( "../helpers/helper.r" )

load( "runs/parameters.Rdata" )

weighted_preds <- as.matrix(read.table("runs/predictions_x.csv",
                              sep = " ",
                              na = "NA"
                              ) )

libs <- read.table("runs/libraries_x.txt",
                   header = FALSE,
                   sep = "\t" )

raw_df <- read.csv(filename,
                   header = TRUE,
                   sep = "," )

truth <- raw_df$x[pred[1]:pred[2]]


mus  <- get_mus(raw_df)
sigs <- get_sigs(raw_df)
df   <- data.frame(scale(raw_df))


