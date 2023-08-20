#' Predicted immuno-prognostic subtypes in single-cell data.
#' @title Predicted immuno-prognostic subtypes in single-cell data.
#' @param data Seurat object.
#' @param split.by Sample class
#' @param method knn, nb, mean, All Default: All (including three models)
#'
#' @return
#' @export
#'
#' @examples
#'test <- predict_scrna(SingleTest,split.by = "orig.ident", method = "All")
predict_scrna <- function(data, #Seurat object
                          split.by = "orig.ident",
                          method = "All" #knn, nb, mean
){
  #
  Tumor_pair <- function(data1){
    hepdata<-data1
    pair2 <- pair_single
    heptest <- lapply(pair2, function(z){
      data <- hepdata
      z1<-data.frame(A=rownames(data)[match(z$C1[,1], rownames(data))],
                     B=rownames(data)[match(z$C1[,2], rownames(data))])
      z1 <-na.omit(z1)
      data <- data[z1[,1],] - data[z1[,2],]
      out<-apply(data, 2, function(x){
        #count <- sum(x==0)
        count <- 0
        sum(x>0)/(length(x)-count)
      })

      data <- hepdata
      z2<-data.frame(A=rownames(data)[match(z$C2[,1], rownames(data))],
                     B=rownames(data)[match(z$C2[,2], rownames(data))])
      z2 <-na.omit(z2)
      data <- data[z2[,1],] - data[z2[,2],]
      out1<-apply(data, 2, function(x){
        #count <- sum(x==0)
        count <- 0
        sum(x>0)/(length(x)-count)
      })
      #
      return(cbind(out,out1))
    })
    heptest<-do.call(cbind,lapply(heptest, function(x) x))
    colnames(heptest) <- c("group1.C1","group1.C2","group2.C1","group2.C2",
                           "group3.C1","group3.C2","group4.C1","group4.C2")
    heptest<- data.frame(heptest)
    return(heptest)
  }

  #
  machine_class <- function(data1,model,split.by = "orig.ident"){
    library(stats)
    count_data <- data1@assays$RNA@counts
    data <- Tumor_pair(count_data)
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
  mean_class<-function(data1,split.by = "orig.ident"){
    count_data <- data1@assays$RNA@counts
    data <- Tumor_pair(count_data)
    data$meanC1 <- apply(data[,c(1,3,5,7)],1,mean)
    data$meanC2 <- apply(data[,c(2,4,6,8)],1,mean)
    cell_subtype <- ifelse(data$meanC1 > data$meanC2, "C1","C2")
    cell <- data.frame(data1[[split.by]],
                       Subtype= cell_subtype)
    ##classification
    sample_name<-as.character(unique(data1[[split.by]])[[split.by]])
    output <- NULL
    temp<-NULL
    for (i in 1:length(sample_name)) {
      temp<-which(data1[[split.by]] == sample_name[i])
      output[i]<- length(which(cell_subtype[temp]=="C1"))   /length(temp)
    }
    sample_ratio <- data.frame(Sample_name=sample_name,
                               Ratio=output,
                               Subtype=ifelse(output>0.5,"C1","C2"))
    out <- list(cell_class=cell,
                sample_class=sample_ratio)
    return(out)
  }

  if(method=="knn"){
    KNN_result <- machine_class(data1=data, model=knn,split.by)
    return(KNN_result)
  }else if(method=="nb"){
    NB_result <- machine_class(data1=data,model=nb, split.by)
    return(NB_result)
  }else if(method=="mean"){
    #
    MEAN_result <- mean_class(data1=data, split.by)
    return(MEAN_result)
  }else{
    #MEAN
    MEAN_result <- mean_class(data1=data, split.by)
    ##KNN
    KNN_result <- machine_class(data1=data, model=knn,split.by)
    ##NB
    NB_result <- machine_class(data1=data,model=nb, split.by)
    Final <-  apply(data.frame(KNN_result$sample_class$Subtype,
                               NB_result$sample_class$Subtype,MEAN_result$sample_class$Subtype), 1, function(x){
      if(sum(x=="C1") >= sum(x=="C2")){
        return("C1")
      }else{
        return("C2")
      }
    })
    out <- list(KNN=KNN_result,NB=NB_result,MEAN=MEAN_result,
                Final_subtype=data.frame(Sample=KNN_result$sample_class$Sample_name,KNN=KNN_result$sample_class$Subtype,
                                         NB=NB_result$sample_class$Subtype,MEAN=MEAN_result$sample_class$Subtype,
                                         Subtype=Final))
    return(out)
  }
}



