#' @title Clustering
#' @description Identifies the cell clusters, i.e. the cell subpopulations.
#'
#' @details If the user knows the number of clusters present in her data set,
#' then `n.cluster` can be set and the estimation of the number of clusters is
#' skipped. `n` is the maximum number of clusters that the automatic estimation
#' of the number of clusters will consider. It is ignored if `n.cluster` is
#' provided. `method` must be "simlr" or "kmeans" exclusively. If set to
#' "simlr", then the function uses the **SIMLR()** function (**SIMLR** package)
#' to perform clustering. If set to "kmeans" the function will perform a
#' dimensionality reduction by principal component analysis (PCA) followed by
#' K-means clustering and 2-dimensional projection by t-distributed stochastic
#' neighbor embedding (t-SNE). Regardless of the value of `method` ("simlr" or
#' "kmeans"), in case `n.cluster` is not provided, then the function relies on
#' the **SIMLR_Estimate_Number_of_Clusters()** function to determine the number
#' of clusters, between 2 and `n`. If `plot` is TRUE, then the function displays
#' the t-SNE map with each cell colored according to the cluster it belongs to.
#' If `method` argument is "simlr", then it further displays a heatmap of the
#' similarity matrix calculated by the **SIMLR()** function. If `pdf` is TRUE,
#' then the function exports the t-SNE plot in a pdf file in the *images*
#' folder. The file is named "t-SNE_map-X.pdf", where X is the `method`
#' argument. If `write` is TRUE, then the function writes two text files in the
#' *data* folder. The first one is called "cluster-Y-X.txt", containing the
#' cluster vector assigning each cell of `data` to a cluster. The second one is
#' called "tsne-Y-X.txt", containing the coordinates of each cell in the 2D
#' t-SNE projection. "X" is the `method` argument anf "Y" is the retained number
#' of clusters.
#' 
#' Note that SIMLR might no longer be available in the most recent versions of R.
#' It is thus necessary to load the library by yourself before calling this function
#' if you want to use it (with \code{library(SIMLR)}).
#'
#' @param data a data frame of n rows (genes) and m columns (cells) of read or UMI counts (note : rownames(data)=genes)
#' @param n.cluster a number, an estimation of the ideal number of clusters is computed if equal to 0
#' @param n a number, the maximum to consider for an automatic determination of the ideal number of clusters
#' @param method "kmeans" or "simlr"
#' @param plot a logical
#' @param pdf a logical
#' @param write a logical
#'
#' @return The function returns a list containing a numeric vector specifying the cluster assignment for each cell,
#' a 2D t-SNE projection, and the number of cells per cluster.
#'
#' @export
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom grDevices dev.off
#' @importFrom grDevices rainbow
#' @importFrom graphics abline
#' @importFrom graphics legend
#' @importFrom graphics symbols
#' @importFrom stats kmeans
#' @importFrom stats prcomp
#' @import Rtsne
#' @import pheatmap
#'
#' @examples
#' data=matrix(runif(100000,0,1),nrow=500,ncol=200)
#' clustering(data,n.cluster=2,method="kmeans")
clustering = function(data,n.cluster=0,n=10,method=c("kmeans","simlr"),plot=TRUE,pdf=TRUE,write=TRUE){
  if (dir.exists("images")==FALSE & pdf==TRUE){
    dir.create("images")
  }
  if (dir.exists("data")==FALSE & write==TRUE){
    dir.create("data")
  }
  if (n.cluster==0){
    cat("Estimating the number of clusters with SIMLR",fill=TRUE)
    if (!requireNamespace("SIMLR", quietly = TRUE)) {
        warning("The SIMLR package must be installed to use this functionality")
        return(NULL)
    }
    c = SIMLR::SIMLR_Estimate_Number_of_Clusters(data, NUMC=2:n)
    a=data.matrix(as.numeric((c$K1+c$K2)))
    rownames(a)=c(2:n)
    n.cluster = as.numeric(rownames(subset(a, a==min(a))))
    cat(paste("Estimated number of clusters =",n.cluster),fill=TRUE)
  }
  method = match.arg(method)

  final=list()
  if (method=="simlr"){
    if (!requireNamespace("SIMLR", quietly = TRUE)) {
      warning("The SIMLR package must be installed to use this functionality")
      return(NULL)
    }
    s = SIMLR::SIMLR(data,c=n.cluster,no.dim = 2)
    Y = s$y
    cluster = Y$cluster
    final[[1]] = cluster
    final[[2]] = s$F
    final[[3]] = Y$size
    final[[4]] = log(s$S*10^6+1)
    names(final) = c("cluster","t-SNE","numbers","similarity")
  }
  if (method=="kmeans"){
    data = data[rowSums(data)>0,]
    pca = prcomp(t(data),center=TRUE,scale.=TRUE)
    v = (pca$sdev^2)/(sum(pca$sdev^2))
    n = min(which(diff(v)>-10^-4))
    if (n==0 | identical(n, integer(0))){
      n=round(ncol(pca$x)/2)
    }
    pca = pca$x[,seq_len(n)]
    tsne = Rtsne(pca)
    km = kmeans(pca,n.cluster)
    cluster = km$cluster
    final[[1]] = cluster
    final[[2]] = tsne$Y
    final[[3]] = km$size
    names(final) = c("cluster","t-SNE","numbers")
  }

  cr=rainbow(max(cluster))
  if (plot==TRUE){
    if (method=="simlr"){
      pheatmap(log(s$S*10^6+1))
    }
    plot(x=final[[2]][,1],y=final[[2]][,2],type='n',main="t-SNE Map",xlab="t-SNE1",ylab="t-SNE2",xlim=c(min(final[[2]][,1])*1.5,max(final[[2]][,1])*1.1))
    abline(h=0)
    abline(v=0)
    symbols(x=final[[2]][,1],y=final[[2]][,2],circles=rep(1,nrow(final[[2]])),inches=0.04,bg=cr[cluster],add=TRUE)
    legend("topleft",legend = paste("cluster",seq_len(max(cluster))),fill = cr,cex = 0.7)
  }
  if (pdf==TRUE){
    pdf(paste("./images/t-SNE_map-",method,".pdf",sep=""))
    plot(x=final[[2]][,1],y=final[[2]][,2],type='n',main="t-SNE Map",xlab="t-SNE1",ylab="t-SNE2",xlim=c(min(final[[2]][,1])*1.5,max(final[[2]][,1])*1.1))
    abline(h=0)
    abline(v=0)
    symbols(x=final[[2]][,1],y=final[[2]][,2],circles=rep(1,nrow(final[[2]])),inches=0.04,bg=cr[cluster],add=TRUE)
    legend("topleft",legend = paste("cluster",seq_len(max(cluster))),fill = cr,cex = 0.7)
    dev.off()
  }
  if (write==TRUE){
    fwrite(data.frame(final[[1]]),paste("./data/cluster-",n.cluster,"-",method,".txt",sep=""),sep="\t")
    fwrite(data.frame(final[[2]]),paste("./data/tsne-",n.cluster,"-",method,".txt",sep=""),sep="\t")
  }
  cat(paste(n.cluster,"clusters detected"),fill=TRUE)
  for (i in seq_len(n.cluster)){
    cat(paste("cluster",i,"->",final[[3]][i],"cells"),fill=TRUE)
  }

  return(final)
}
