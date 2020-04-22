#' @title Data Prepare
#' @description Prepares the data for further analysis
#'
#' @details `file` is the path to the file containing the read or UMI count
#' matrix the user wants to analyze.
#' @details
#' `most.variables` can be set to N to select the Nth most variables genes. This
#' option allows the user to use a reduced matrix (N x number of cells) to
#' perform the clustering step faster.
#' @details
#' `lower` and `upper` are used to remove the genes whose average counts are outliers.
#' The values of these arguments are fractions of the total number of genes and
#' hence must be between 0 and 1. Namely, if `lower = 0.05`, then the function
#' removes the 5% less expressed genes and if `upper = 0.05`, then the function
#' removes the 5% most expressed genes.
#' @details
#' If `normalize` is FALSE, then the function skips the 99th percentile
#' normalization and the log transformation.
#' @details
#' If `write` is TRUE, then the function writes two text files. One for the
#' normalized and gene thresholded read counts table and another one for the
#' genes that passed the lower and upper threshold. Note that the length of the
#' genes vector written in the *genes.txt* file is equal to the number of rows
#' of the table of read counts written in the *data.txt* file.
#'
#' @param file a string for the scRNAseq data file
#' @param most.variables a number
#' @param lower a number in [0,1], low quantile threshold
#' @param upper a number in [0,1], high quantile threshold
#' @param normalize a logical, if TRUE, then computes 99th percentile normalization
#' @param write a logical
#' @param verbose a logical
#' @param plot a logical
#'
#' @return The function returns a data frame of filtered and/or normalized data
#' with genes as row names.
#'
#' @export
#'
#' @import data.table
#' @importFrom stats quantile
#' @import graphics
#' @importFrom stats sd quantile
#'
#'
#' @examples
#' file <- system.file("scRNAseq_dataset.txt",package = "SingleCellSignalR")
#' data <- data_prepare(file = file)
data_prepare <- function(file, most.variables=0, lower=0, upper=0,normalize=TRUE,write=FALSE,verbose=TRUE, plot=FALSE){

  if (dir.exists("data")==FALSE & write==TRUE){
    dir.create("data")
  }
  if (!file.exists(file)){
    stop(paste(file,"doesn't exist."))
  }

  data <- fread(file,data.table=FALSE)

  for (i in seq_len(ncol(data))){
    if (is.character(data[,i])==FALSE){
      p=i-1
      break
    }
  }
  if (sum(data[,p] %in% c(mm2Hs$`Mouse gene name`,mm2Hs$`Gene name`))==0){
    stop("Please convert gene ID's (Ensembl or NCBI ID) to official HUGO
        gene symbols")
  }
  data <- subset(data, !duplicated(data[,p]))
  genes <- data[,p]
  data <- data[,-c(seq_len(p))]
  rownames(data) <- genes
  data <- data[rowSums(data)>0,]
  data <- data.frame(data[,apply(data,2,function(x) quantile(x,0.99))>0])
  if (normalize==TRUE){
    cat("log-Normalization",fill=TRUE)
    q <- apply(data,2,quantile,0.99)
    data <- log(1+sweep(data,2,q/median(q),"/"))
  }
  data <- data[rowSums(data)>0,]
  data <- data[rowSums(data)<=quantile(rowSums(data),1-upper) &
                rowSums(data)>=quantile(rowSums(data),lower),]

  if (verbose==TRUE){
    cat(paste(dim(data)[1],"genes"),fill=TRUE)
    cat(paste(dim(data)[2],"cells"),fill=TRUE)
    cat(paste("Zero rate = ",round(sum(data==0)*1000/prod(dim(data)))/10,"%",sep=""),fill=TRUE)
  }
  if (write==TRUE){
    fwrite(data.frame(data),"./data/data.txt",sep="\t")
    fwrite(data.frame(rownames(data)),"./data/genes.txt",sep="\t")
  }
  if (most.variables!=0){
    m <- apply(data,1,mean)
    cv <- apply(data,1,sd)/m
    names(cv) <- rownames(data)
    cv <- cv[m>quantile(m,0.5)]
    if (length(cv)<most.variables){
      mv.genes <- names(cv)
    } else {
      mv.genes <- names(sort(cv,decreasing=TRUE))[seq_len(most.variables)]
    }
    mv.genes <- names(sort(cv,decreasing=TRUE))[seq_len(most.variables)]
    res <- list(data,data[mv.genes,])
    names(res) <- c("complete.dataset","most.var.dataset")
  } else {
    res <- data
  }
  return(res)
}

