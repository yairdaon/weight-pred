#!/usr/bin/Rscript
library( rEDM )
source( "../helpers/plotting.r" )
df <- read.csv("original.csv")
df$time <- NULL
df <- scale(df)
noise <- rnorm( prod(dim(df)), mean=0, sd=sqrt(0.1))  
df    <- noise + df ## As long as first_column_time == FALSE

for( i in 1:ncol(df) )
{
    print(i)
    ts <- df[ ,i]
    output <- simplex(ts,
                      E = 1:10,
                      lib = c(1000,3000),
                      stats_only = TRUE)
    x11()
    plot(output$E,
         output$rho,
         type = "l")
    hold()
}

    
