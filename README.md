# IPS29GPS
The goal of IPS29GPS is to predict immuno-prognostic subtypes

## Installation

You can install the development version of IPS29GPS from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Sixiegg/IPS29GPS")
```

## Example

This is a basic example:

``` r
#library(IPS29GPS)
## Predicted immuno-prognostic subtypes in bulk transcriptomic data.
bulkresult <- predict_29pairs(bulktest,Pairs,genetype = "Symbol")
head(bulkresult)

library(IPS29GPS)
library(GSVA)
library(Seurat)
## Predicted immuno-prognostic subtypes in single-cell data.

# ssGSEA result
ssGSEA_result <- predict_scrna(SingleTest,gene_set,split.by="orig.ident",method="ssGSEA")
str(ssGSEA_result)

##RF result
RF_result <- predict_scrna(SingleTest,gene_set,split.by="orig.ident",method="RF")
str(RF)

##KNN result
KNN_result <- predict_scrna(SingleTest,gene_set,split.by="orig.ident",method="KNN")
str(KNN)

##Final result
Final_result <- predict_scrna(SingleTest,gene_set,split.by="orig.ident",method="All")
str(Final_result)
#final subtype
Final_result$Final_subtype

```
