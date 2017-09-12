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


show_ccm <- function(var,
                     tp,
                     libsizes=1:10 * 20,
                     n_samples=100,
                     E=4 )
{
    
    var_x_chl_rhos <- var_xmap_chl(var,
                                   tp,
                                   libsizes,
                                   n_samples,
                                   E )
    
    chl_x_var_rhos <- chl_xmap_var(var,
                                   tp,
                                   libsizes,
                                   n_samples,
                                   E )

    tit <- paste0( "Cross maps of Chlorophyll and ", var )

    x11()
    plot(libsizes,
         var_x_chl_rhos,
         main = tit,
         type = "l",
         col = "blue",
         ylim = c(-1,1) )

    lines(chl_x_var_rhos,
          col = "red")
         
    hold(1)
    

    ## print( paste0( "Chl X ", var, " tp = ", tp, ", rho = ", chl_xmap_var( var, tp = -2 )$rho ) )
    ## print( paste0(  var, " X Chl, tp = ", tp, ", rho = ", var_xmap_chl( var, tp = -2 )$rho ) )
}

