#!/usr/bin/Rscript
## title: "Huisman"
## author: "Ethan with George & Yair"
## date: "11/20/2017"
## output: html_document
## ---

## ```{r setup, include=FALSE}
## knitr::opts_chunk$set(echo = TRUE)
## ```


## packages
## library(tidyverse)
library(deSolve)
## library(gridExtra)
## library(gridGraphics)
## library(dplyr)

## First we define a generating function for the model. Note that a
## second set of parameters is included that reproduces another version
## of the system from the same text.

## model generator
gen_Huisman <- function(n, tau = 10)
{
    ## initial values & coefficients (Huisman & Weissing 2001, Fig. 2)
    N0 <- c(N1=0.1, N2=0.1, N3=0.1, N4=0.1, N5=0.1)
    R0 <- c(R1=10, R2=10, R3=10)
    
    num_N <- length(N0)
    num_R <- length(R0)
    
    
    ps <- list()
    ps$S <- R0
    ps$K <- matrix(c(0.20, 0.05, 0.50, 0.05, 0.50,
                     0.15, 0.06, 0.05, 0.50, 0.30,
                     0.15, 0.50, 0.30, 0.06, 0.05), nrow = num_R, byrow = TRUE)
    ps$C <- matrix(c(0.20, 0.10, 0.10, 0.10, 0.10,
                     0.10, 0.20, 0.10, 0.10, 0.20,
                     0.10, 0.10, 0.20, 0.20, 0.10), nrow = num_R, byrow = TRUE)
    ##   ps$K <- matrix(c(0.20, 0.05, 1.00, 0.05, 1.20, 
    ##                 0.25, 0.10, 0.05, 1.00, 0.40, 
    ##                 0.15, 0.95, 0.35, 0.10, 0.05), nrow = num_R, byrow = TRUE)
    ##   ps$C <- matrix(c(0.20, 0.10, 0.10, 0.10, 0.10, 
    ##                 0.10, 0.20, 0.10, 0.10, 0.20, 
    ##                 0.10, 0.10, 0.20, 0.20, 0.10), nrow = num_R, byrow = TRUE)# step function
    ps$r <- rep.int(1, num_N) ## maximum growth rate
    ps$m <- rep.int(0.25, num_N) ## mortality
    ps$D <- 0.25
    

    
    dF <- function(t,x,ps)
    {
        N <- x[1:num_N]
        R <- x[num_N+(1:num_R)]
        mu <- ps$r * apply(R / (ps$K + R), 2, min)
        
        return(list(c(N * (mu - ps$m) ,
                      ps$D * (ps$S - R) - colSums(t(ps$C) * mu * N))))
    }

    timepoints <- seq(1,n)*tau
    ##model_data <- lsoda(c(N0,R0),timepoints,dF,parms =  ps,hmax = 0.01,maxsteps = 5000)
    model_data <- ode(c(N0,R0),timepoints,dF,parms =  ps,hini = 0.01,"rk4")
    model_data <- as.data.frame(model_data)
    
    ## return(list(N = N, R = R))
    return(model_data[,1:(num_N+1)])
    ## return(model_data)
}

## Generate a realization of the model.
block_huisman <- gen_Huisman(3000)
write.csv(block_huisman,
          "originals/huisman.csv",
          quote = FALSE,
          row.names = FALSE )








######################################
## Rest ##############################
######################################

## ## block_huisman <- mutate( block_huisman , regime=as.numeric(N2 > N4))

## ## To see visualize how the different variables "see" and distort the
## ## native attractor, we define 3 different immersions of the data.

## ## define embeddings
## L_views <- list(
##     list(vars=c("N1","N2","N4"),delays=c(0,0,0)),
##     list(vars=c("N2","N2","N2"),delays=0:-2),
##     list(vars=c("N4","N4","N4"),delays=0:-2)
## )

## ## Source the function "make_block" copied from past projects.


## make_block <- function(data, cols, delays, 
##                        lib = c(1, NROW(data)), 
##                        diff_col = rep(FALSE, length(cols)))
## {   ## Takes an input matirx or data frame and creates a block of lag coordinates
##     ## to use with rEDM.
##     ## INPUTS:
##     ##   data - matrix, array, or data.frame with time series variables arranged in columns
##     ##   cols - vector indices or names of the columns of 'data' to use for each column of block
##     ##   delays - vector with same length as cols specifying the time displacement for that column
##     ##   diff_col - vector of logical on whether to apply first differencing to this lag coordinate variable
##     ##
##     ## OUTPUT:
##     ##   block - array with length(cols) of columns, and NROW(data) rows
    
##     if(!is.numeric(cols)){
##         cols <- as.numeric(factor(levels=colnames(data),x = cols))
##     }
    
##     lib <- matrix(lib,ncol = 2)
##     ## data <- as.matrix(data)
    
##     ncol <- length(cols)
##     nrow <- dim(data)[1]
##     block <- as.data.frame(array(NA,dim = c(nrow,ncol)))
##     names(block) <- 1:ncol
    
##     for (i in 1:ncol)
##     {
##         I <- 1:nrow
##         I_delay <- intersect(I,I+delays[i])
##         block[I_delay-delays[i],i] <- data[I_delay,cols[i]]
##         if (delays[i] < 0){
##             ## remove data points that fall at start of lib segments
##             block[lib[,1] - (0:(delays[i]+1)),i] <- NA
##             names(block)[i] <- paste(colnames(data)[cols[i]],'_t-',abs(delays[i]),sep="")  
##         } else if (delays[i] > 0) {
##             ## remove data points that fall at end of lib segments
##             block[lib[,2] - (0:(delays[i]-1)),i] <- NA
##             names(block)[i] <- paste(colnames(data)[cols[i]],'_t+',abs(delays[i]),sep="")  
##         } else {
##             names(block)[i] <- paste(colnames(data)[cols[i]],'_t',sep="")
##         }
        
##         if (diff_col[i]){
##             block[,i] <- c(NA,diff(block[,i]))
##         }
##     } ## i
    
    
##     return(block)
## }

## ## We define a function that makes a 3D plot panel.

## plot_attractor <- function(vars,delays){
##     block_temp <- make_block(data = block_huisman,cols = vars,delays=delays)
##     pdf("huisman.pdf")
##     plot3D::lines3D(block_temp[,1],
##                     block_temp[,2],
##                     block_temp[,3],
##                     colvar=block_huisman$regime)
    
##     dev.off()
##     ## grid.echo()
##     ## g <- grid.grab()
##     ## return(g)
## }

## ## print( L_views[[1]] )
## ## plot_attractor( L_views[[1]]$vars, L_views[[1]]$delays )

## ## ## Finally make a set of 3 plots.
## ## lapply(L_views,function(view){ 
## ##     plot_attractor(view$vars,view$delays)
## ## }),
## ## nrow=1)

## ## r generate 3 panel

## ## do.call(grid.arrange,
## ##         c(
## ##             lapply(L_views,function(view){ 
## ##                 plot_attractor(view$vars,view$delays)
## ##             }),
## ##             nrow=1)
## ##         )
