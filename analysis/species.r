#!/usr/bin/env Rscript
library(rEDM)

df <- read.table("originals/data_20111121.txt",
                  header = TRUE,
                  sep = "\t",
                  ## nrows = 100,
                  na.strings ="NaN" )
dates <- as.POSIXlt(as.Date(paste(df$Year, df$Month, df$Day, sep="-")))


## Inferring how much Chlorophyll every cell is worth

## A linear model of Chl-A against cell counts
sum_lm <- lm( df$Chlorophyll..surface ~
                  df$All.Diatoms..cell.counts +
                  df$All.Dinoflagellates..cell.counts )

## Geerating the inferred Chl levels from the model.
cof <- coef( sum_lm )
dino_diatom <-
    cof[1] +
    cof[2] * df$All.Diatoms..cell.counts +
    cof[3] * df$All.Dinoflagellates..cell.counts


## Plotting:

## Save to this file
pdf( "plots/species.pdf", width = 50 )

## Which data to show
ind <- !is.na( dino_diatom )

## Chlorophyll abundance
plot(dates[ind],
     df$Chlorophyll..surface[ind],
     col = "black",
     xlab = "Time",
     ylab = "Abundance",
     type = "p",
     main = "Measured Chl and inferred Chl from diatom and dinoflagellate measurements"
     )

## Diatom + Dinoflagellates data.
points(dates[ind],
       dino_diatom[ind],
       col = "red"
       )

## 95% Threshold.
abline( quantile(df$Chlorophyll..surface, 0.95, na.rm = TRUE ), 0, col = "green" )

## Legends... Boring...
legend("topleft",
       legend = c( "Chl", "Diatoms + Dinoflagellates", "Chl 95% Threshold" ),
       lty = c(1,1,1),
       col = c( "black", "red", "green" ))
legend("topright",
       legend = c( "Chl", "Diatoms + Dinoflagellates", "Chl 95% Threshold" ),
       lty = c(1,1,1),
       col = c( "black", "red", "green" ))
legend("center",
       legend = c( "Chl", "Diatoms + Dinoflagellates", "Chl 95% Threshold" ),
       lty = c(1,1,1),
       col = c( "black", "red", "green" ))

## Finalize saving to file.
dev.off()





## Junk code

## all_lm <- lm(
##     df$Chlorophyll..surface ~
##         df$Alexandrium.spp..cell.counts +
##         df$Cochlodinium.spp..cell.counts +
##         df$Dinophysis.spp..cell.counts +
##         df$Gymnodinium.spp..cell.counts +
##         df$Lingulodinium.polyedrum.cell.counts +
##         df$Prorocentrum.spp..cell.counts +
##         df$Pseudo.nitzschia.spp..cell.counts +
##         ## df$Non.toxic.Diatoms..cell.counts +
##         ## df$Non.toxic.Dinoflagellates..cell.counts +
##         ## df$Total.non.toxic.cells..cell.counts +
##         ## df$Total.cells..5um..cell.counts +
##         ## df$All.Diatoms..cell.counts +
##         ## df$All.Dinoflagellates..cell.counts +
##         df$Prorocentrum.micans.cell.counts +
##         df$Akashiwo.sanguinea.cell.counts +
##         df$Pseudo.nitzschia.seriata.group.cell.counts +
##         df$Pseudo.nitzschia.delicatissima.group.cell.counts
## )


## cof <- coef( all_lm )
## all_n_all <-
##     cof[1] +
##     cof[2] * df$Alexandrium.spp..cell.counts +
##     cof[3] * df$Cochlodinium.spp..cell.counts +
##     cof[4] * df$Dinophysis.spp..cell.counts +
##     cof[5] * df$Gymnodinium.spp..cell.counts +
##     cof[6] * df$Lingulodinium.polyedrum.cell.counts +
##     cof[7] * df$Prorocentrum.spp..cell.counts +
##     cof[8] * df$Pseudo.nitzschia.spp..cell.counts +
##     cof[9] * df$Prorocentrum.micans.cell.counts +
##     cof[10] * df$Akashiwo.sanguinea.cell.counts +
##     cof[11] * df$Pseudo.nitzschia.seriata.group.cell.counts +
##     cof[12] * df$Pseudo.nitzschia.delicatissima.group.cell.counts
