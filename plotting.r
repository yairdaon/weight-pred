hold <- function(n=1){
    message("Press Return To Continue")
    invisible(readLines("stdin", n))
}

weight_plot <- function(time,
                        obs,
                        weights,
                        filename,
                        xtension )
    
{
    jpeg( paste0( "plots/", filename, "_", xtension, ".jpg" ) )

    par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
    plot(time,
         obs,
         type = "l",
         col  = "black",
         xlab = "time(days)",
         ylab = "chlorophyll A",
         main = paste0( "Chlorophyll observations and ", gsub( "_", " ", filename), " weights.")
         )
    par(new = TRUE)
    plot(time,
         weights,
         type = "l",
         col = "red",
         axes = FALSE,
         bty = "n",
         xlab = "",
         ylab = "")
    axis( side=4, at = pretty(range(weights)) )
    mtext("weights", side=4, line=3 )
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

