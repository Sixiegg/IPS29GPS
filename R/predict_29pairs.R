#' Predicted immuno-prognostic subtypes in bulk transcriptomic data.
#'
#' @param data Composed of gene symbol or id and exp.
#' @param Pairs 29 gene pairs
#' @param genetype Symbol or ID
#'
#' @return
#' @export
#'
#' @examples
#'predict_29pairs(bulktest,Pairs,genetype = "Symbol")
#'
predict_29pairs <- function(data,Pairs,
                            genetype = "Symbol" ##"ID")
){
  #
  geneid <- data[,1]
  data <- data[,-1]
  Aexp <- data[match(Pairs[[genetype]][,1],geneid),]
  Bexp <- data[match(Pairs[[genetype]][,2],geneid),]
  Cexp <- na.omit(Bexp - Aexp)
  print(paste0("A total of ",dim(Cexp)[1]," pairs were retained!!!"))
  Count <- colSums(Cexp>=0)
  Subtype <- ifelse(Count >= dim(Cexp)[1]/2, "C1", "C2")
  Ratio <- Count / dim(Cexp)[1]
  out <- data.frame(Sample = colnames(Cexp), Subtype=Subtype, Ratio=Ratio)
  return(out)
}
