#!/usr/bin/Rscript
source( "../helpers/plotting.r" )


neighbors <- 1:4
cols <- c( "red", "blue", "black", "green" )


for( method in c( "mve", "uwe" ) )
{
    load( paste0( "runs/", method, "_parameters.Rdata" ) )
    pdf(paste0("plots/", method, ".pdf" ) )
    plot(0,
         0,         
         main = paste0( toupper(method), "Skill for Different Numbers of Nearest Neighbors"),
         xlab = "Library size",
         ylab = "Skill(rho)",
         xlim = c(min(lib_sizes),max(lib_sizes)),
         ylim = c(0.5,1) 
         )
    
    for( i in neighbors )
    {
        df <- read.csv(paste0("runs/", method, "_y_", i, ".csv"), header = TRUE, sep = "," )
        
        lines(df$lib_sizes,
              df$mean,
              type = "l",
              lty = 1,
              col = i)
        lines(df$lib_sizes,
              df$bot,
              lty = 3,
              col = i)
        lines(df$lib_sizes,
              df$top,
              lty = 3,
              col = i)
        
    }
    
    grid(lwd = 1, lty = 1)
    
    legend(x = "bottomright",
           title = "# nn",
           legend = neighbors,
           col = cols,
           lwd = 1)
    
    dev.off()

}
