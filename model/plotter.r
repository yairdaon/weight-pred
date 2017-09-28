#!/usr/bin/Rscript
source( "../helpers/plotting.r" )
df <- read.csv("data/rhos.csv", header = TRUE, sep = "," )


pdf("plots/weighted_prediction.pdf")
plot(df$lib_sizes,
     df$avg,
     type = "l",
     lty = 1,
     col = "red",
     ylim = c(0.6,1),
     main = "Uncertainty-Weighted Prediction Skill for y with Quartiles",
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
