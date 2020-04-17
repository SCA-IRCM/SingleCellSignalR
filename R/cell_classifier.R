#' @title Cell classifier
#' @description  Classifies cells using cell type specific markers.
#'
#' @details
#' The ` markers` argument must be a table with cell type gene signatures, one
#' cell type in each column. The column names are the names of the cell types.
#' @details
#' The *markers.default* table provides an example of this format.
#' @details
#' If ` tsne` is not provided, then the function will just not display the cells
#' on the t-SNE. Although t-SNE maps are widely used to display cells on a 2D
#' projection, the user can provide any table with two columns and a number of
#' rows equal to the number of columns of `data` (e.g. the two first components
#' of a PCA).
#' @details
#' If ` plot.details` is TRUE, then the function plots the number of cells
#' attributed to a single cell type as a function of the threshold applied to
#' the normalized gene signature average.
#' @details
#' If ` write` is TRUE, then the function writes four different text files. (1)
#' The "raw classification matrix" provides the normalized average gene
#' signature for each cell type in each individual cell, a number between 0 and
#' 1. This matrix has one row per cell type and one column per cell, and the
#' sum per column is 1. Row names are the cell type names (column names of the
#' markers table) and the column names are the individual cell identifiers
#' (column names of `data`). (2) The "thresholded classification matrix", which
#' is obtained by eliminating all the values of the "raw classification matrix"
#' that are below a threshold a\*. In practice, a\* is automatically determined
#' by the function to maximize the number of cells that are assigned to a single
#' cell type and all the cells (columns) assigned to 0 or >1 cell types are
#' discarded. The number of cells assigned to a single type depending on a\*
#' can be plotted by using the parameter `plot.details=TRUE`. (3) A cluster
#' vector assigning each cell to a cell type. Note that a supplementary,
#' virtual cluster is created to collect all the cells assigned to 0 or >1
#' types. This virtual cluster is named "undefined". (4) A table associating
#' each cell type to a cluster number in the cluster vector.
#'
#' @param data a data frame of n rows (genes) and m columns (cells) of read or
#' UMI counts (note : rownames(data)=genes)
#' @param genes a character vector of HUGO official gene symbols of length n
#' @param markers a data frame of cell type signature genes
#' @param tsne (optional) a table of n rows and 2 columns with t-SNE
#' projection coordinates for each cell
#' @param plot.details a logical (if TRUE, then plots the number of cells
#' attributed to one cell type, see below)
#' @param write a logical
#' @param verbose a logical
#'
#' @return The function returns a list containing the thresholded table, the
#' maximum table, the raw table, a cluster vector and the cluster names. The
#' maximum table is a special thresholded table where in every column only the
#' maximum gene signature is kept. It can be used to force the classification
#' of every cell.
#'
#' @export
#' @import pheatmap
#' @importFrom graphics plot
#' @importFrom data.table melt
#'
#' @examples
#' data <- matrix(runif(1000,0,1),nrow=50,ncol=20)
#' rownames(data) <- paste("gene",seq_len(50))
#' markers <- matrix(paste("gene",seq_len(10)),ncol=5,nrow=2)
#' colnames(markers)=paste("type",seq_len(5))
#' cell_classifier(data,rownames(data),markers)
cell_classifier <- function(data,genes,markers=markers_default,tsne=NULL,
                           plot.details=FALSE,write=TRUE,verbose=TRUE){
  if (dir.exists("cell-classification")==FALSE & write==TRUE){
      dir.create("cell-classification")
  }
  rownames(data) <- genes
  n.types <- ncol(markers)

  tmp <- data[as.character(unlist(markers))[as.character(unlist(markers)) %in%
                                             genes],]

  if (is.null(dim(tmp))==TRUE){
    cat("Not enough markers genes to pursue the cell classification",fill=TRUE)
    return()
  }
  final <- matrix(0,ncol=ncol(tmp),nrow=(ncol(markers)))

  for (i in seq_len(n.types)){
    m.genes <- markers[,i][markers[,i] %in% genes]
    final[i,] <- colSums(tmp[m.genes,])/length(m.genes)
  }
  final[is.na(final)] <- 0
  final[,colSums(final)!=0] <- apply(final[,colSums(final)!=0],2,function(x)
    x/sum(x))
  # final <- round(final*1000)/10
  rownames(final) <- colnames(markers)
  colnames(final) <- colnames(data)
  l <- matrix(0,101,2)
  q <- 0
  for (n in seq(0.01,1,0.01)){
    q <- +1
    f <- final
    f[f<n] <- 0
    l[q+1,1] <- n
    l[q+1,2] <- sum(apply(f,2,function(x) sum(x==0)==(n.types-1)))
  }
  seuil <- max(l[which(l[,2]==max(l[,2])),1])
  m <- matrix(0,(n.types+1),2)
  for (n in 0:n.types){
    f <- final
    f[f<seuil] <- 0
    f <- matrix(f[,apply(f,2,function(x) sum(x==0)==n)],nrow=n.types)
    m[n+1,1] <- n
    m[n+1,2] <- sum(apply(f,2,function(x) sum(x==0)==n))
  }
  m <- matrix(m[m[,2]!=0,],ncol=2)

  res <- list()
  final.s <- final
  final.s[final.s<seuil] <- 0
  res[[1]] <- final.s[,apply(final.s,2,function(x) sum(x==0)==(n.types-1))]
  res[[3]] <- final
  final[!apply(final,2,function(x) x==max(x))] <- 0
  res[[2]] <- final
  h <- pheatmap::pheatmap(res[[1]],cluster_rows=TRUE,cluster_cols=TRUE,
                         show_colnames=FALSE)

  nn <- seq_len(length(rownames(res[[1]])[rowSums(res[[1]])!=0]))
  names(nn) <- rownames(res[[1]])[rowSums(res[[1]])!=0]
  cluster <- unlist(lapply(apply(final.s,2,function(x) names(x[x>0])),
                          function(x) if(length(x)==1){
                            return(nn[as.character(x)])} else {return(0)}))
  cluster[cluster==0]=max(cluster)+1

  if (sum(m[,2]!=0)==1){
    n <- names(nn)
  } else {
    n <- c(names(nn), "Undefined cells")
  }
  cr <- c(rainbow(max(cluster)-1),"gray")
  res[[4]] <- cluster
  # res[[5]] <- h
  res[[5]] <- c(rownames(res[[1]][rowSums(res[[1]])>0,]),"undefined")
  names(res) <- c("tresh_mat","max_mat","raw_mat","cluster","c.names")
  d <- data.frame(Cell_type=c(rownames(res[[1]]),"undefined"),
               Number_of_cells=c(apply(res[[1]],1,function(x) sum(x!=0)),
                                 ncol(res[[3]])-ncol(res[[1]])))
  d <- d[d$Number_of_cells!=0,]
  rownames(d) <- paste("cluster",seq_len(nrow(d)))


  if (verbose==TRUE){
    cat(paste("A threshold of",round(seuil*1000)/10,"% maximizes the number of
              cells assigned to one cell type"),fill=TRUE)
    for (i in (nrow(m)):1){
      cat(paste(m[i,2],"cell(s) or",round(m[i,2]*1000/ncol(data))/10,
                "% identified to",(n.types)-m[i,1],"cell type(s)"),sep=" ",
          fill=TRUE)
    }
    cat(" ",fill=TRUE)
    print(d)
  }

  if (is.null(tsne)==FALSE){
    plot(x=tsne[,1],y=tsne[,2],type='n',main="t-SNE Map",xlab="t-SNE1",
         ylab="t-SNE2",xlim=c(min(tsne[,1])*2,max(tsne[,1])*1.05))
    abline(h=0)
    abline(v=0)
    symbols(x=tsne[,1],y=tsne[,2],circles=rep(1,nrow(tsne)),inches=0.04,
            bg=cr[cluster],add=TRUE)
    legend("topleft",legend=n,fill=cr,cex=0.75)
  }
  if (plot.details==TRUE){
    plot(l,xlab="Threshold (%)",ylab="Number of cells",
         main= "Attribution of one cell to one cell type",type='l')
    abline(v=l[l[,1]==seuil,1], col="red", lty=2)
    symbols(l[l[,1]==seuil,1],l[l[,1]==seuil,2],circles=1,inches=0.05,
            bg="red",fg="red",add=TRUE)
    text(x=l[l[,1]==seuil,1],y=l[l[,1]==seuil,2]-0.07*(max(l)),
         labels=paste(seuil,"%"),offset=20)
  }
  if (write==TRUE){
    fwrite(data.frame(cbind(rownames(res[[1]]),res[[1]])),
           "cell-classification/threshold_class_matrix.txt",sep="\t")
    fwrite(data.frame(cbind(rownames(res[[3]]),res[[3]])),
           "cell-classification/raw_class_matrix.txt",sep="\t")
    nn <- cbind(names(nn),nn)
    colnames(nn) <- c("cluster.name","cluster.number")
    fwrite(data.frame(nn),"cell-classification/cluster2name.txt",sep="\t")
    fwrite(data.frame(cluster),"cell-classification/cluster.txt",sep="\t")
    c.names <- res[[5]]
  }

  return(res)
}
