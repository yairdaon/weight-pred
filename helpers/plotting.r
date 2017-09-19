hold <- function(n=1){
    message("Press Return To Continue")
    invisible(readLines("stdin", n))
}


double_plot <- function(x,  ## time
                        y1, ## chlorophyll
                        y2, ## uncertainty
                        q = quantile( y1, 0.95, na.rm = TRUE),
                        xlabel   = "x",
                        y1label  = "black",
                        y2label  = "red",
                        title    = "",
                        filename = FALSE
                        )
    
{
    if( filename == FALSE )
        x11()
    else
        pdf( paste0("plots/", filename, ".pdf") )
    
    par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for y2 axis
    plot(x,
         y1,
         type = "l",
         col  = "black",
         xlab = xlabel,
         ylab = y1label,
         main = title
         )
    abline(q, 0, col = "green" )
        
    par(new = TRUE)
    plot(x,
         y2,
         type = "l",
         col = "red",
         axes = FALSE,
         bty = "n",
         xlab = "",
         ylab = "")
    axis( side=4, at = pretty(range(y2)) )
    mtext(y2label, side=4, line=3 )

    if( filename == FALSE )
        hold(1)
    else
        dev.off()
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

