## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(out.width = "100%")
knitr::opts_chunk$set(fig.align = 'center')
library(knitr)
library(dataSDA)
library(RSDA)

## -----------------------------------------------------------------------------
data(mushroom)
head(mushroom)

## -----------------------------------------------------------------------------
mushroom.set <- set_variable_format(data = mushroom, location = 8, var = "Species")
head(mushroom.set)

## -----------------------------------------------------------------------------
mushroom.tmp <- RSDA_format(data = mushroom.set, sym_type1 = c("I", "S"),
                            location = c(25, 31), sym_type2 = c("S", "I", "I"),
                            var = c("Species", "Stipe.Length_min", "Stipe.Thickness_min"))
head(mushroom.tmp)

## -----------------------------------------------------------------------------
mushroom.clean <- clean_colnames(data = mushroom.tmp)
head(mushroom.clean)

## ---- eval = FALSE------------------------------------------------------------
#  write_csv_table(data = mushroom.clean, file = 'mushroom_interval.csv')

## -----------------------------------------------------------------------------
mushroom.int <- read.sym.table(file = 'mushroom_interval.csv', header = T, sep = ';', dec = '.', row.names = 1)
head(mushroom.int)

## -----------------------------------------------------------------------------
data(Abalone.iGAP)
head(Abalone.iGAP)

## -----------------------------------------------------------------------------
Abalone <- iGAP_to_MM(Abalone.iGAP, c(1, 2, 3, 4, 5, 6, 7))
head(Abalone)

