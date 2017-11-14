#!/usr/bin/Rscript
library(rEDM)
source( "helpers/plotting.r" )

## Consider all possible sums of to models, check if it predicts
## better.

args <- commandArgs(trailingOnly = TRUE)
args <- paste( unlist(args), collapse = " " )
if( length(args) != 0 ) {

    if( grepl("oos", args, ignore.case = TRUE ) ) 
        xtension <- "OoS"
    if( grepl("loocv", args, ignore.case = TRUE ) ) 
        xtension <- "LOOCV" 

    load(paste0("data/", tolower(xtension), "_full_run.Rdata"))
    
    if( grepl("rho", args, ignore.case = TRUE ) )
    {
        metric <- function( ts1, ts2 ) 
            return( cor(ts1, ts2, use = "complete.obs" ) )
        xtension <- paste0( xtension, "_rho" )
    } else if( grepl("mae", args, ignore.case = TRUE ) )
    {
        metric <- function( ts1, ts2 ) 
            return(-mean(abs(ts1-ts2), na.rm = TRUE    ) )
        xtension <- paste0( xtension, "_mae" )
    } 

} else {
    stop( "No command line arguments. Aborting." )
}

dims     <- dim( pred_table )
n_models <- dims[1]
n        <- dims[2]

metrics <- numeric( n_models )
for( i in 1:n_models )
    metrics[i] <- metric(chl_p1wk, pred_table[i, ])


perms  <- combn( n_models, 2 )
n_comb <- dim( perms )[2]

## We always maximize, so we start at negative infinity
best <- -Inf 


for( i in 1:n_comb)
{

    if(!(i %% 100000))
        print(i)

    ## Get model indices
    model1 <- perms[1,i]
    model2 <- perms[2,i]
    
    ## Calculate what their sum prediction is
    sum_preds <- pred_table[model1, ] + pred_table[model2, ]
    
    ## Find  metric value for the summed prediction
    val <- metric(chl_p1wk,sum_preds)
                  
    ## If it is better then previous best - update.
    if( val > best )
    {
        models <- c(model1, model2)
        best <- val
    }

}


print( paste0( "Metric of summed predictor = ", best,  "." ) )

desc1 <- paste( other_vars[ combinations[ , models[1] ] ], collapse = ", " ) 
pred1 <- pred_table[ models[1], ]
val1  <- metric(chl_p1wk, pred1)
print( desc1 )
print( paste0( "Metric = ", val1,  "." ) )

desc2 <- paste( other_vars[ combinations[ , models[2] ] ], collapse = ", " )
pred2 <- pred_table[ models[2], ]
val2  <- metric(chl_p1wk, pred2)
print( desc2 )
print( paste0( "Metric = ", val2,  "." ) )


svg( paste0("plots/additive_", tolower(xtension), ".svg"), width = 100 ) 
plot(chl_p1wk,
     type='l', 
     bty = "l",
     col = "black",
     ylim = c(-1, 10)
     )
lines(pred1 + pred2,
      col = "red")
lines(pred1,
      col = "blue")
lines(pred2,
      col = "green")
abline(norm_threshold,
       0,
       col = "yellow" )
legend("topleft",
       legend = c( "True Chl", "Summed predictors", desc1, desc2, "95% threshold" ),
       lty = 1,
       bty = 'n',
       cex = .75,
       col = c( "black", "red", "blue", "green", "yellow" ) )
dev.off()

print("There's nothing to do, there's nowhere to go, I need action, I need adventure, I need to punch!!!")
