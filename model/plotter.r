#!/usr/bin/Rscript
source( "../helpers/plotting.r" )


neighbors <- 1
cols <- c( "red", "blue", "black", "green" )

pdf(paste0("plots/preds.pdf" ) )

plot(0,
     0,         
     main = "Skill for Different Methods",
     xlab = "Library size",
     ylab = "Skill(rho)",
     xlim = c(10,100),
     ylim = c(0.5,1) 
     )

df <- read.csv(paste0("runs/hao_y_", neighbors, ".csv"),
               header = TRUE,
               sep = "," )
lines(df$lib_sizes,
      df$mean,
      type = "l",
              lty = 1,
      col = cols[1])


df <- read.csv(paste0("runs/mve_y_", neighbors, ".csv"),
               header = TRUE,
               sep = "," )
lines(df$lib_sizes,
      df$mean,
      type = "l",
      lty = 1,
      col = cols[2])

df <- read.csv(paste0("runs/uwe_y_", neighbors, ".csv"),
               header = TRUE,
               sep = "," )
lines(df$lib_sizes,
      df$mean,
      type = "l",
      lty = 1,
      col = cols[3])



grid(lwd = 1, lty = 1)

legend(x = "bottomright",
       title = "Method",
       legend = c( "Hao", "mine", "UWE" ),
       col = cols[1:3],
       lwd = 1)

dev.off()




## lines(df$lib_sizes,
##       df$bot,
##       lty = 3,
##       col = i)
## lines(df$lib_sizes,
##       df$top,
##       lty = 3,
##       col = i)
        
    
