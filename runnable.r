#!/usr/bin/Rscript
library(rEDM)
library(zoo)
source( "../helpers/helper.r" )

## Get a weighted prediction predictors and their uncertainty
mve_pred <- function(pred_table,
                     ranks)
{
    ord <- order( rhos, decreasing = TRUE )
    ind <- ord[ 1 : ceiling(sqrt(length(ranks))) ]
    ret <- colMeans( pred_table[ ind, ] )
    return( ret )
}

## Get a uncertainty-weighted prediction
uwe_pred <- function(pred_table, var_table, method = "exp")
{
    ## Get rid of those zero variance predictions.
    var_table[ var_table == 0 ] <- NA
    
    ## Get increasing order, so lower uncertainties
    ## have lower indices.
    ord <- colOrder( var_table )

    ## Take best sqrt(k), according to the MVE heuristics
    best <- ceiling(sqrt(nrow(var_table)))
    ind  <- ord[1:best, ]

    if( method == "exp" )
        w <- exp(-matrix(var_table [ind], nrow = best, ncol = ncol(var_table) ) )
    else if( method == "minvar" )
        w <- 1 / matrix(var_table [ind], nrow = best, ncol = ncol(var_table) )
    else
        stop( paste0("Unknown weighting scheme '", method "'." ) )
    
    p <- matrix(pred_table[ind], nrow = best, ncol = ncol(var_table) )

    return( colSums( w*p ) / colSums(w) )
    
}

dummy_df <- data.frame( lib_sizes = numeric(0),
	    		mean = numeric(0),
			bot  = numeric(0),
			med  = numeric(0),
 			top  = numeric(0))
write.table(dummy_df,
    	    file = "data/uwe.csv",
            sep = ",",
            append = FALSE,
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
write.table(dummy_df,
    	    file = "data/mve.csv",
            sep = ",",
            append = FALSE,
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)



## If you want to download the file for some reason... don't!
## filename = paste0("http://science.sciencemag.org/highwire/filestream/683325/",
##                   "field_highwire_adjunct_files/1/aag0863_SupportingFile_Suppl._Excel_seq1_v2.xlsx")

## Load data and immediately normalize to zero mean and unit std dev
df <-normalize_df(read.csv("originals/three_species.csv",
                           header = TRUE,
                           sep = "," )
                  )

## Parameters
n_vars <- ncol(df) ## == 3 cuz x, y, z
E <- 3
n_lags <- 3 ## Lags: 0, -1, -2
n_samp <- 100 ## Number of random libraries
pred <- c(2501, 3000) ##c(2501,3000) ## Prediciton always constant
pred_size <- pred[2] - pred[1] + 1
lib_sizes <- (1:20)*5

## Craeate lags for every variable and make y_p1
df <- lag_every_variable(df, n_lags)
y_p1 <- c( as.vector(lag(zoo(df$y),1)) , NA )
df <- data.frame(y_p1 = y_p1, df )

## We know the allowed combinations must have one unlagged
## variable. By the way function combn is structured, these will be
## the first n_comb (defined below).
combinations <- combn( n_vars*n_lags, E ) ## Get ALL cobminations.
## Only this many are OK.
n_comb <- choose( n_vars*n_lags, E )- choose(n_vars*(n_lags-1), E ) 
## Throw away the rest.
combinations <- combinations[ ,1:n_comb ]

## Preallocate memory for predictions and uncertainties
pred_table <- matrix( NA,  nrow = n_comb, ncol = pred_size )
var_table  <- matrix( Inf, nrow = n_comb, ncol = pred_size )
rhos       <- matrix( NA,  nrow = n_samp, ncol = 2         )
prediction <- numeric(                           pred_size )
ranks      <- numeric(            n_comb                   )

for (lib_size in lib_sizes )
{
    ## For every random library starting point
    for( smp in 1:n_samp )
    {
        
        ## Choose random lib of length exactly 100. Pred is always the same
        ## and they never overlap.
        lib  <- sample(501:2000, 1)
        lib  <- c(lib, lib + lib_size - 1 )
	if( smp %% 5 == 0 )
 	       print( paste0("Sample ", smp, ", library size == ", lib_size ) )
        
	for( i in 1:n_comb )
        {
            cols <- combinations[,i]
            
            ## Just to make sure, all cols have one unlagged coordinate.
            stopifnot( min(cols) <= n_vars )
            cols <- cols + 1 ## cuz first col is y_p1
            
            output <- block_lnlp(df,
                                 lib = lib,   ## lib is chosen randomly above
                                 pred = pred, ## pred is ALWAYS the same
                                 method = "simplex",
                                 tp = 0, ## Zero step cuz it is built into y_p1
                                 columns = cols, 
                                 target_column = "y_p1",
                                 first_column_time = FALSE,
                                 stats_only = FALSE )
            
            pred_table[i, ] <- output$model_output[[1]]$pred
            var_table [i, ] <- output$model_output[[1]]$pred_var

            ranks[i] <- block_lnlp(df,
                                   lib = lib,   ## lib is chosen randomly above
                                   pred = lib, ## pred is ALWAYS the same
                                   method = "simplex",
                                   tp = 0, ## Zero step cuz it is built into y_p1
                                   columns = cols, 
                                   target_column = "y_p1",
                                   first_column_time = FALSE,
                                   stats_only = TRUE )$rho
       
        }
        
        truth <- output$model_output[[1]]$obs
        uwe <- uwe_pred( pred_table, var_table)
        mve <- mve_pred( pred_table, ranks ) 
        rhos[smp,1] <- cor(uwe, truth, use = "complete.obs" )
        rhos[smp,2] <- cor(mve, truth, use = "complete.obs" )
    }

    probs <- c(0.25,0.5,0.75)
    uwe_vec <- matrix( c( lib_size, mean(rhos[,1]), quantile(rhos[,1], probs = probs) ), nrow = 1)
    mve_vec <- matrix( c( lib_size, mean(rhos[,2]), quantile(rhos[,2], probs = probs) ), nrow = 1)
    
    write.table(uwe_vec,
    	        file = "data/uwe.csv",
                sep = ",",
                append = TRUE, 
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
    write.table(mve_vec,
    	        file = "data/mve.csv",
                sep = ",",
                append = TRUE, 
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
    
}

uwe_df <- read.csv("data/uwe.csv", header = TRUE, sep = "," )
mve_df <- read.csv("data/mve.csv", header = TRUE, sep = "," )

pdf("plots/predictions.pdf")
plot(uwe_df$lib_sizes,
     uwe_df$avg,
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
      mve_df$avg,
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

legend(x = "topleft",
       legend = c("UWE", "MVE" ),
       col = c( "red", "blue" ),
       lwd = 1)
dev.off()
