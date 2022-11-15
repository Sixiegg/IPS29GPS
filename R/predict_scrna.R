#' Predicted immuno-prognostic subtypes in single-cell data.
#' @title Predicted immuno-prognostic subtypes in single-cell data.
#' @param data Seurat object. Please use FindVariableFeatures function from Seurat to identify highly variable features !!!
#' @param gene_set A list of C1 and C2 gene sets
#' @param split.by Sample class
#' @param method ssGSEA, RF, KNN. Default: All (including three models)
#'
#' @return
#' @export
#'
#' @examples
#'
predict_scrna <- function(data, #Seurat object
                          gene_set,
                          split.by = "orig.ident",
                          method = "All" #ssGSEA, RF, KNN
){
  ssGSEA_class <- function(data,
                           result,
                           gene_set,
                           split.by = "orig.ident"){
    scale_row<-function (x)
    {
      m = apply(x, 1, mean, na.rm = T)
      s = apply(x, 1, sd, na.rm = T)
      return((x - m)/s)
    }
    result1<-scale_row(result)
    C2 <- result1[2,,drop=F]
    C1 <- result1[1,,drop=F]
    result2<-(C1-C2)>0
    sample_name<-as.character(unique(data[[split.by]])[[split.by]])
    ##classification
    output <- NULL
    temp<-NULL
    for (i in 1:length(sample_name)) {
      temp<-which(data[[split.by]] == sample_name[i])
      output[i]<-sum(result2[temp])/length(temp)
    }
    cell <- data.frame(data[[split.by]],
                       C1_scale_score= as.numeric(C1) ,
                       C2_scale_score= as.numeric(C2) ,
                       Subtype= c(ifelse(result2,"C1","C2") )
    )
    sample_ratio <- data.frame(Sample_name=sample_name,
                               Ratio=output,
                               Subtype=ifelse(output>0.5,"C1","C2"))
    out <- list(cell_class=cell,
                sample_class=sample_ratio)
    return(out)
  }
  #
  machine_class <- function(data1,model,split.by = "orig.ident"){
    pair_model<-do.call(rbind,lapply(strsplit(model[["coefnames"]],"_"), function(x) x))
    data<-data1@assays$RNA@counts
    data <- data[pair_model[,1],]-data[pair_model[,2],]
    data[]<-apply(data,2,function(x) ifelse(x>0,1,x))
    data[]<-apply(data,2,function(x) ifelse(x<0,-1,x))
    rownames(data)<-model[["coefnames"]]
    data<-t(as.matrix(data))
    library(stats)
    cell_subtype <- stats::predict(model, newdata = data)

    sample_name<-as.character(unique(data1[[split.by]])[[split.by]])
    ##classification
    output <- NULL
    temp<-NULL
    for (i in 1:length(sample_name)) {
      temp<-which(data1[[split.by]] == sample_name[i])
      output[i]<- length(which(cell_subtype[temp]=="C1"))   /length(temp)
    }
    cell <- data.frame(data1[[split.by]],
                       Subtype= cell_subtype )
    sample_ratio <- data.frame(Sample_name=sample_name,
                               Ratio=output,
                               Subtype=ifelse(output>0.5,"C1","C2"))
    out <- list(cell_class=cell,
                sample_class=sample_ratio)
    return(out)
  }

  count_data <- data@assays$RNA@counts
  if(length(data@assays$RNA@var.features)==0){
    message('Please use FindVariableFeatures function from Seurat to identify highly variable features !!!')
    return(NA)
  }

  if(method=="ssGSEA"){
    count_data <- count_data[match(data@assays$RNA@var.features,rownames(count_data)),]
    count_score <- gsva(count_data,gene_set,method='ssgsea',kcdf="Poisson")
    ssGSEA_result <- ssGSEA_class(data,count_score,gene_set,split.by)
    return(ssGSEA_result)
  }else if(method=="RF"){
    RF_result <- machine_class(data=data, model=RF, split.by)
    return(RF_result)
  }else if(method=="KNN"){
    KNN_result <- machine_class(data=data, model=KNN, split.by)
    return(KNN_result)
  }else{
    ##ssGSEA
    count_data <- count_data[match(data@assays$RNA@var.features,rownames(count_data)),]
    count_score <- gsva(count_data,gene_set,method='ssgsea',kcdf="Poisson")
    ssGSEA_result <- ssGSEA_class(data,count_score,gene_set,split.by)
    ##RF
    RF_result <- machine_class(data=data, model=RF, split.by)
    ##KNN
    KNN_result <- machine_class(data=data, model=KNN, split.by)
    Final <-  apply(data.frame(ssGSEA_result$sample_class$Subtype,
                               RF_result$sample_class$Subtype,KNN_result$sample_class$Subtype), 1, function(x){
      if(sum(x=="C1") >= sum(x=="C2")){
        return("C1")
      }else{
        return("C2")
      }
    })
    out <- list(ssGSEA=ssGSEA_result,RF=RF_result,KNN=KNN_result,
                Final_subtype=data.frame(Sample=RF_result$sample_class$Sample_name,
                                         Subtype=Final))
    return(out)
  }
}



