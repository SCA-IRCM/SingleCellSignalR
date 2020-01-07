#' @title most variable interactions
#' @description Displays a heatmap showing the most variable interactions over
#' all clusters.
#'
#'
#' @param data a data frame of n rows (genes) and m columns (cells) of read or
#' UMI counts (note : rownames(data)=genes)
#' @param genes a character vector of HUGO official gene symbols of length n
#' @param cluster a numeric vector of length m
#' @param c.names (optional) cluster names
#' @param n an integer the number of most variables interactions
#' @param species "homo sapiens" or "mus musculus"
#'
#' @return The function displays a heatmap showing the most variable
#' interactions over all clusters
#' @export
#'
#' @importFrom pheatmap pheatmap
#' @importFrom stats var
#'
#' @examples
#' data=matrix(runif(1000,0,1),nrow=5,ncol=200)
#' genes = c("gene 1","gene 2","gene 3","gene 4","gene 5")
#' cluster=c(rep(1,100),rep(2,100))
#' mv_interactions(data,genes,cluster)
mv_interactions = function(data,genes,cluster,c.names=NULL,n=30,
                           species=c("homo sapiens","mus musculus")){
  if (is.null(c.names)==TRUE){
    c.names = paste("cluster",seq_len(max(cluster)))
  }
  species = match.arg(species)
  if (species=='mus musculus'){
    Hs2mm = mm2Hs[,1]
    mm2Hs = mm2Hs[,2]
    names(mm2Hs) = Hs2mm
    names(Hs2mm) = as.character(mm2Hs)
    m.names = mm2Hs[genes]
    data = subset(data,(!is.na(m.names)))
    m.names = m.names[!is.na(m.names)]
    genes=as.character(m.names)
  }
  rownames(data) = genes
  l_sc = LRdb[is.element(LRdb$ligand,rownames(data)),]
  int_sc = l_sc[is.element(l_sc$receptor,rownames(data)),]
  lr_sc=matrix(0,nrow=nrow(int_sc),ncol=max(cluster)^2)
  rownames(lr_sc) = paste(int_sc$ligand,int_sc$receptor,sep=" / ")
  med = sum(data)/(nrow(data)*ncol(data))
  nam=NULL
  q=0
  for (i in seq_len(max(cluster))){
    for (j in seq_len(max(cluster))){
      q=q+1
      lr_sc[,q] = (rowMeans(data[int_sc$ligand,cluster==i])*
                     rowMeans(data[int_sc$receptor,cluster==j]))^0.5/
        (med + (rowMeans(data[int_sc$ligand,cluster==i])*rowMeans(
          data[int_sc$receptor,cluster==j]))^0.5)
      nam=c(nam,paste(c.names[i],c.names[j],sep=" -> "))
    }
  }
  colnames(lr_sc) = nam
  if (sum(lr_sc)!=0){
    if (nrow(lr_sc)<n){
      n = nrow(lr_sc)
    }
    lr_sc=subset(lr_sc,rowSums(lr_sc)!=0)
    v=apply(lr_sc,1,var)/apply(lr_sc,1,mean)
    lr_sc = lr_sc[order(v,decreasing = TRUE),]
    lr_sc = lr_sc[apply(lr_sc,1, max)>0.5,]
    pheatmap::pheatmap(lr_sc[seq_len(n),colSums(lr_sc[seq_len(n),])!=0],
                       cluster_cols = TRUE)
  } else {
    cat("No interactions detected. Make sure the genes vector is composed of
        HUGO official gene names.",fill=TRUE)
  }

}


