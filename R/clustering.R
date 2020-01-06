#' @title Clustering
#' @description Identifies the cell clusters, i.e. the cell subpopulations.
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
#' @import SIMLR
#' @import Rtsne
#' @import pheatmap
#'
#' @examples
#' data=matrix(runif(100000,0,1),nrow=500,ncol=200)
#' clustering(data,n.cluster=2,method="kmeans")
clustering = function(data,n.cluster=0,n=10,method=c("simlr","kmeans"),plot=TRUE,pdf=TRUE,write=TRUE){
  if (dir.exists("images")==FALSE & pdf==TRUE){
    dir.create("images")
  }
  if (dir.exists("data")==FALSE & write==TRUE){
    dir.create("data")
  }
  if (n.cluster==0){
    cat("Estimating the number of clusters",fill=TRUE)
    c = SIMLR_Estimate_Number_of_Clusters(data, NUMC=2:n)
    a=data.matrix(as.numeric((c$K1+c$K2)))
    rownames(a)=c(2:n)
    n.cluster = as.numeric(rownames(subset(a, a==min(a))))
    cat(paste("Estimated number of clusters =",n.cluster),fill=TRUE)
  }
  method = match.arg(method)

  final=list()
  if (method=="simlr"){
    s = SIMLR(data,c=n.cluster,no.dim = 2)
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
