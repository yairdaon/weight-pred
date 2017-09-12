#!/usr/bin/Rscript
library(rEDM)
source( "helpers/plotting.r" )

args <- commandArgs(trailingOnly = TRUE)
if( (length(args) == 0)  ) {
    xtension <- "Test"
} else if (grepl("oos", args, ignore.case = TRUE )) {
    xtension <- "OoS"
} else {
    xtension <- "LOOCV"
}
load(paste0("data/", tolower(xtension), "_full_run.Rdata"))
dims <- dim( pred_table )
n <- dims[2]

var_table[ is.na(pred_table) ] <- NA
avg_var <- colMeans( var_table , na.rm = TRUE )

double_plot(1:n,
            chl,
            avg_var,
            q = norm_threshold,
            xlabel   = "time",
            y1label  = "ChlA",
            y2label  = "uncertainty - variance",
            title    = paste0("Prediction uncertainty and Chl-A with 95% Chl-A threshold, ", xtension), 
            filename = paste0( "var_chl_", tolower(xtension) )
            )
