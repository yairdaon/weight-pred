#!/usr/bin/Rscript
source( "../helpers/plotting.r" )

uwe_df <- read.csv("runs/uwe_y.csv", header = TRUE, sep = "," )
mve_df <- read.csv("runs/mve_y.csv", header = TRUE, sep = "," )

pdf("plots/predictions.pdf")
plot(uwe_df$lib_sizes,
     uwe_df$mean,
     type = "l",
     lty = 1,
     col = "red",
     ylim = c(0.2,1),
     main = "Skill for UWE and MVE",
     xlab = "Library size",
     ylab = "skill(rho)")
lines(uwe_df$lib_sizes,
      uwe_df$bot,
      lty = 3,
      col = "red")
lines(uwe_df$lib_sizes,
      uwe_df$top,
      lty = 3,
      col = "red")


lines(mve_df$lib_sizes,
      mve_df$mean,
      lty = 1,
      col = "blue")
lines(mve_df$lib_sizes,
      mve_df$bot,
      lty = 3,
      col = "blue")
lines(mve_df$lib_sizes,
      mve_df$top,
      lty = 3,
      col = "blue")

grid(lwd = 1, lty = 1)

legend(x = "topleft",
       legend = c("UWE", "MVE" ),
       col = c( "red", "blue" ),
       lwd = 1)

dev.off()
