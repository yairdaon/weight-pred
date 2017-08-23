library( zoo )

hold <- function(n=1){
    message("Press Return To Continue")
    invisible(readLines("stdin", n=1))
}

harry_plotter <- function(x, tit)
{
    tit <- paste( "Transformed ", tit, " data" )
    x <- na.omit( x )
    x11()
    hist(x,
         main = tit,
         breaks = 20,
         probability=TRUE,
         xlim = c(-4, 4 ) )
    
    lines( density(x), col="red" )
    z <- seq(-4,4, length = 1000)
    lines(z,
          dnorm(z, mean=0, sd = 1 ),
          col = "blue" )
    
    hold()
}
