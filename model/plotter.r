#!/usr/bin/Rscript
source( "../helpers/plotting.r" )

var_names <- c( "x", "y", "z" )
for( var_name in var_names )
{

    df <- read.csv(paste0("data/rhos_", var_name, ".csv"), header = TRUE, sep = "," )

    pdf(paste0( "plots/weighted_prediction_", var_name, ".pdf") )
    plot(df$lib_sizes,
         df$avg,
         type = "l",
         lty = 1,
         col = "red",
         ylim = c(0.6,1),
         main = paste0("Uncertainty-Weighted Prediction Skill for ",var_name, " with Quartiles"),
         xlab = "Library size",
         ylab = "skill(rho)")
    lines(df$lib_sizes,
          df$bot,
          lty = 3,
               col = "red")
    lines(df$lib_sizes,
          df$top,
          lty = 3,
          col = "red")
    dev.off()
}
