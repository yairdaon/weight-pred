#!/usr/bin/Rscript
library(rEDM)
library(zoo)
source( "../helpers/helper.r" )

## Get a weighted prediction predictors and their uncertainty
mve_pred <- function(pred_table,
                rhos)
{
    ord <- order( rhos, decreasing = TRUE )
    ind <- ord[ 1 : ceiling(sqrt(length(rhos))) ]
    ret <- colMeans( pred_table[ ind, ] )
    return( ret )
}

## Get a uncertainty-weighted prediction
uwe_pred <- function(pred_table, var_table)
{
    ## Get rid of those zero variance predictions.
    var_table[ var_table == 0 ] <- NA
    
    ## Get increasing order, so lower uncertainties
    ## have lower indices.
    ord <- colOrder( var_table )

    ## Take best sqrt(k), according to the MVE heuristics
    best <- ceiling(0.15*(nrow(var_table)))
    ind <- ord[1:best, ]
    
    w <- exp(-matrix(var_table [ind], nrow = best, ncol = ncol(var_table) ) )
    ## w <- 1  / matrix(var_table [ind], nrow = best, ncol = ncol(var_table) )
    p <-      matrix(pred_table[ind], nrow = best, ncol = ncol(var_table) )

    return( colSums( w*p ) / colSums(w) )
    
}


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
n_samp <- 6 ##100 ## Number of random libraries
pred <- c(2501, 3000) ##c(2501,3000) ## Prediciton always constant
pred_size <- pred[2] - pred[1] + 1
lib_sizes <- (3:5)*10

## Library changes in size and is also random
uwe_results <- data.frame(lib_sizes = lib_sizes)
mve_results <- data.frame(lib_sizes = lib_sizes)

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
rhos       <- matrix( NA,  nrow = 2,      ncol = n_samp    )
prediction <- numeric(                           pred_size )
ranks      <- numeric(            n_comb                   )

lib_ind <- 1
for (lib_size in lib_sizes )
{
    ## For every random library starting point
    for( smp in 1:n_samp )
    {
      
        
        ## Choose random lib of length exactly 100. Pred is always the same
        ## and they never overlap.
        lib  <- sample(501:2000, 1)
        lib  <- c(lib, lib + lib_size - 1 )
        print( paste0("Sample ", smp, ", library ", lib[1], "-", lib[2] ) )
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
            
            val <- output$model_output[[1]]$pred
            err <- output$model_output[[1]]$pred_var

            stopifnot( length(val) == pred[2] -pred[1] + 1 )
            stopifnot( length(err) == pred[2] -pred[1] + 1 )
            pred_table[i, ] <- val
            var_table [i, ] <- err

            ranks[i] <- block_lnlp(df,
                                   lib = lib,   ## lib is chosen randomly above
                                   pred = lib, ## pred is ALWAYS the same
                                   method = "simplex",
                                   tp = 0, ## Zero step cuz it is built into y_p1
                                   columns = cols, 
                                   target_column = "y_p1",
                                   first_column_time = FALSE,
                                   ##short_output = TRUE,
                                   stats_only = TRUE )$rho
       
        }
        
        truth <- output$model_output[[1]]$obs
        uwe <- uwe_pred( pred_table, var_table)
        mve <- mve_pred( pred_table, ranks ) 
        rhos[1, smp] <- cor(uwe, truth, use = "complete.obs" )
        rhos[2, smp] <- cor(mve, truth, use = "complete.obs" )
    }

    quarts <- quantile(rhos[1, ], probs = c(0.25,0.5,0.75))
    uwe_results$avg[lib_ind] <- mean(rhos[1, ])
    uwe_results$bot[lib_ind] <- quarts[1]
    uwe_results$med[lib_ind] <- quarts[2]
    uwe_results$top[lib_ind] <- quarts[3]

    quarts <- quantile(rhos[2, ], probs = c(0.25,0.5,0.75))
    mve_results$avg[lib_ind] <- mean(rhos[2, ])
    mve_results$bot[lib_ind] <- quarts[1]
    mve_results$med[lib_ind] <- quarts[2]
    mve_results$top[lib_ind] <- quarts[3]

    lib_ind <- lib_ind + 1

    
}

write.csv(uwe_results,
          file = "data/uwe.csv",
          row.names = FALSE,
          quote = FALSE,
          na = "NA")

write.csv(mve_results,
          file = "data/mve.csv",
          row.names = FALSE,
          quote = FALSE,
          na = "NA" )

uwe_df <- read.csv("data/uwe.csv", header = TRUE, sep = "," )
mve_df <- read.csv("data/mve.csv", header = TRUE, sep = "," )

print( "UWE, then MVE:" )
print( uwe_df$avg )
print( mve_df$avg )

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
