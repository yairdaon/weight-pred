hold <- function(n=1){
    message("Press Return To Continue")
    invisible(readLines("stdin", n))
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

