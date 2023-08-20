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
library(Seurat)
## Predicted immuno-prognostic subtypes in single-cell data.
##KNN result
KNN_result <- predict_scrna(SingleTest, split.by="orig.ident",method="knn")

##NB result
NB <- predict_scrna(SingleTest, split.by="orig.ident",method="nb")


##MEAN result
MEAN_result <- predict_scrna(SingleTest, split.by="orig.ident",method="mean")
str(MEAN)

##Final result
Final_result <- predict_scrna(SingleTest, split.by="orig.ident",method="All")

#final subtype
Final_result$Subtype

```
