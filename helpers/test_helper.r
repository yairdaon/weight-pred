#!/usr/bin/Rscript
source( "helper.r" )
##############################################
## Tests for function descale
df <- data.frame(x = c(1,4,6,2,4,6,-1),
                 y = c(1,2,3,4,5,6,7),
                 z = c(0,1,0,1,0,1,0))
df1 <- scale( df )
df2 <- descale(df1)
stopifnot(all(
    abs(df2-df) < 1e-14
))

df3 <- descale(df1,
               mu  = attr(df1, "scaled:center"),
               sig = attr(df1, "scaled:scale" )  
               )
stopifnot(all(
    abs(df3-df) < 1e-14
))
print( "Finished test descale" )


###################################################
## date2serial and serial2date tests, if u wanna use them
ex_date <- "2008-02-29"
ex_day <- date2serial( ex_date ) 
rec_date <- serial2date( ex_day )
rec_day <- date2serial( rec_date )
stopifnot( ex_date == rec_date )
stopifnot( ex_day  == rec_day  )
print( "Finished test serial2date date2serial" )

#####################################################
## Test name_combinations
df <- data.frame(x = numeric(1),
                 y = numeric(1),
                 z = numeric(1))
test_comb <- name_combinations(df, 2, 3)
## print(name_comb)



######################################################
## Test colOrder
m <- 6
n <- 9
X <- matrix(sample(1:20, m*n, replace = TRUE ), nrow = m, ncol = n )
X[ 4 ] <- NA 
X[ 9 ] <- NA 
Y <- matrix(X[colOrder(X)], ncol = ncol(X))

## Check Y columns are sorted and that they have same values in
## corresponding X columns
for( i in 1:ncol(Y) )
{
    stopifnot( !is.unsorted(Y[,i], na.rm = TRUE ) )
    stopifnot( all (Y[,i] %in% X[,i]) )
    stopifnot( all (X[,i] %in% Y[,i]) )
}
print( "Finished test colOrder" )

#######################################################
## Test respect lib
## df <- data.frame(x = c(1,4,5,8,7,8,4,2,5,2,5,6 ),
##                  y = c(5,7,3,9,3,2,5,1,0,8,5,6 ),
##                  z = c(2,5,2,6,2,5,9,7,4,7,2,8 ))
## max_lag <- 4
## n_vars <- ncol(df)
## lib <- c(6,8)

## df <- lag_every_variable(df, max_lag)
## print(df)
## df <- respect_lib( df, lib, max_lag )
## print(df)



comparator <- function( df1, df2 )
{
    df1[ is.na(df1) ] <- -1
    df2[ is.na(df2) ] <- -1
    return( all.equal(df1, df2) )
}


##########################################################
## Tests for lag_every_variable
df <- data.frame(x = c(1,4,5,8,7,8,4,2,5,2,5,7 ),
                 y = c(5,7,3,9,3,2,5,1,0,8,4,6 ))
lib <- matrix(c(5,7),  ncol = 2, byrow = TRUE)
max_lag <- 2

lag_one_test <- data.frame(
    ##TIME = c( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11,12),
    x      = c( 1, 4, 5, 8, 7, 8, 4, 2, 5, 2, 5, 7 ),
    y      = c(5,  7, 3, 9, 3, 2, 5, 1, 0, 8, 4, 6 ),
    x_1    = c(NA, 1, 4, 5,NA, 7, 8, 4, 2, 5, 2, 5 ),
    y_1    = c(NA, 5, 7, 3,NA, 3, 2, 5, 1, 0, 8, 4 )
)
lag_one <- lag_every_variable(df, max_lag = 2)
lag_one <- respect_lib(lag_one, lib, max_lag)
stopifnot( comparator(lag_one,lag_one_test) )

print( "Finished test lag_every_variable and respect_lib" )
## lag_two_test <- data.frame(
##     TIME   = c( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11,12),
##     x      = c( 1, 4, 5, 8, 7, 8, 4, 2, 5, 2, 5, 7 ),
##     y      = c( 5, 7, 3, 9, 3, 2, 5, 1, 0, 8, 4, 6 ),
##     x_1    = c(NA,NA, 1, 4, NA,NA,7, 8, 4, 2, 5, 2 ),
##     y_1    = c(NA,NA, 5, 7, NA,NA,3, 2, 5, 1, 0, 8 ) )
## lag_two <- make_block(df, max_lag = 2, t = NULL, lib = lib, tau = 2)
## stopifnot( comparator( lag_two, lag_two_test ) )


## lag_three_test <- data.frame(
##     TIME   = c( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12 ),
##     x      = c( 1, 4, 5, 8, 7, 8, 4, 2, 5, 2, 5, 7 ),
##     y      = c( 5, 7, 3, 9, 3, 2, 5, 1, 0, 8, 4, 6 ),
##     x_1    = c(NA,NA,NA, 1,NA,NA,NA, 7, 8, 4, 2, 5 ),
##     y_1    = c(NA,NA,NA, 5,NA,NA,NA, 3, 2, 5, 1, 0 )
## )
## lag_three <- make_block(df, max_lag = 2, t = NULL, lib = lib, tau = 3)
## stopifnot( comparator( lag_three, lag_three_test ) )
