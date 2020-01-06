#' @title Data Prepare
#' @description Prepares the data for further analysis
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
#' @return The function returns a data frame of filtered and/or normalized data with genes as row names. If most.variables is equal to N (>0)
#' the function returns also a sub matrix of the N most variables genes.
#' @export
#'
#' @import data.table
#' @importFrom stats quantile
#' @import graphics
#' @importFrom stats sd quantile
#'
#'
#' @examples
#' \donttest{data_prepare("~/scRNAseq_dataset.txt")}
data_prepare = function(file, most.variables=0, lower=0, upper=0,normalize=TRUE,write=TRUE,verbose=TRUE, plot=FALSE){

  if (dir.exists("data")==FALSE & write==TRUE){
    dir.create("data")
  }

  data = fread(file,data.table = FALSE)

  for (i in seq_len(ncol(data))){
    if (is.character(data[,i])==FALSE){
      p=i-1
      break
    }
  }
  if (sum(data[,p] %in% c(mm2Hs$`Mouse gene name`,mm2Hs$`Gene name`))==0){
    cat("Please convert gene ID's (Ensembl or NCBI ID) to official HUGO
        gene symbols",fill=TRUE)
    return()
  }
  data = subset(data, !duplicated(data[,p]))
  genes = data[,p]
  data = data[,-c(seq_len(p))]
  rownames(data) = genes
  data = data[rowSums(data)>0,]
  data = data.frame(data[,apply(data,2,function(x) quantile(x,0.99))>0])
  if (normalize==TRUE){
    cat("log-Normalization",fill=TRUE)
    q = apply(data,2,quantile,0.99)
    data = log(1+sweep(data,2,q/median(q),"/"))
  }
  data = data[rowSums(data)>0,]
  data = data[rowSums(data)<quantile(rowSums(data),1-upper) &
                rowSums(data)>quantile(rowSums(data),lower),]

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
    m = apply(data,1,mean)
    cv = apply(data,1,sd)/m
    names(cv) = rownames(data)
    cv = cv[m>quantile(m,0.5)]
    if (length(cv)<most.variables){
      mv.genes = names(cv)
    } else {
      mv.genes = names(sort(cv,decreasing = TRUE))[seq_len(most.variables)]
    }
    mv.genes = names(sort(cv,decreasing = TRUE))[seq_len(most.variables)]
    res = list(data,data[mv.genes,])
    names(res) = c("complete.dataset","most.var.dataset")
  } else {
    res = data
  }
  return(res)
}

