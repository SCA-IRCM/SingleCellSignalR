#' @title Cluster Analysis
#' @description Analysis of the differentially expressed genes in the clusters
#' and their composition by a marker based approach.
#'
#' @param data a data frame of n rows (genes) and m columns (cells) of read or
#' UMI counts (note : rownames(data)=genes)
#' @param genes a character vector of HUGO official gene symbols of length n
#' @param cluster a numeric vector of length m
#' @param c.names a vector of cluster names
#' @param dif.exp a logical (if TRUE, then computes the diferential gene
#' expression between the clusters using **edgeR**)
#' @param s.pval a value, a fixed p-value threshold
#' @param markers a table of cell type signature genes
#' @param write a logical
#' @param verbose a logical
#'
#' @return The function returns a list comprised of a table of differentially
#' expressed genes, a table of cell types, and a
#' table of cell cluster types.
#'
#' @export
#' @import edgeR
#' @importFrom limma makeContrasts
#' @importFrom stats model.matrix
#'
#' @examples
#' data=matrix(runif(1000,0,1),nrow=5,ncol=200)
#' rownames(data) = c("A2M","LRP1","AANAT","MTNR1A","ACE")
#' cluster=c(rep(1,100),rep(2,100))
#' cluster_analysis(data,rownames(data),cluster,dif.exp=FALSE)
cluster_analysis = function(data,genes,cluster,c.names=NULL,dif.exp=TRUE,
                            s.pval=10^-2,markers=NULL,write=TRUE,verbose=TRUE){
  if (dir.exists("cluster-analysis")==FALSE){
    dir.create("cluster-analysis")
  }
  if (is.null(c.names)==TRUE){
    c.names = paste("cluster",seq_len(max(cluster)))
  }
  if (min(cluster)!=1){
    cluster = cluster + 1 - min(cluster)
  }
  if (length(c.names)!=max(cluster) | sum(duplicated(c.names))>0 |
      grepl("/",paste(c.names,collapse =""))){
    cat("The length of c.names must be equal to the number of
        clusters and must contain no duplicates. The cluster names must not
        include special characters")
    return()
  }
  rownames(data) = genes
  n.cluster = max(cluster)
  z=c(seq_len(n.cluster))
  class = list()
  length(class) = n.cluster
  diff.genes=list()
  final.2=NULL
  types=NULL
  k=0
  v=NULL
  if (dif.exp==TRUE){
    dge = DGEList(data,genes=genes)
    dge = calcNormFactors(dge)
    if (verbose==TRUE){
      cat("edgeR differential gene expression (dge) processing:",fill=TRUE)
    }

    for (i in z){
      if (verbose==TRUE){
        cat(paste("Looking for differentially expressed genes in", c.names[i]),
            fill=TRUE)
      }
      cl = as.matrix(cluster)
      cl[cluster==i,] = 1
      cl[cluster!=i,] = 2
      c = factor(cl)
      design = model.matrix(~0+c)
      cm = makeContrasts(c1-c2,levels=design)
      y = estimateDisp(dge,design,robust=TRUE)
      fit.y = glmFit(y,design)
      lrt = glmLRT(fit.y,contrast=cm)
      sel.r = topTags(lrt,p.value=s.pval,adjust.method="BH",n=nrow(y$counts))
      if (length(sel.r)==0){
        if (verbose==TRUE){
          cat(paste("No differentially expressed genes in", c.names[i]),
              fill=TRUE)
        }
        v = c(v,i)
      } else {
        k = k+1
        resu = sel.r$table
        resu = resu[order(resu$logFC,decreasing = TRUE),]
        diff.genes[[k]] = resu
        if (write==TRUE){
          fwrite(data.frame(resu),paste("./cluster-analysis/table_dge_",
                                        c.names[i],".txt",sep=""),sep="\t")
        }
      }
    }
    if (is.null(v)==TRUE){
      names(diff.genes) = c.names
    } else {
      names(diff.genes) = c.names[-v]
    }
  }

  if (is.null(markers)==FALSE){
    final = matrix(0,nrow=n.cluster,ncol=ncol(markers))
    colnames(final) = colnames(markers)
    rownames(final) = c.names
    for (i in z){
      tmp=NULL
      for (j in seq_len(ncol(markers))){
        m.genes = markers[,j][markers[,j] %in% genes]
        tmp = c(tmp,sum(rowSums(data.frame(data[m.genes,
                                                cluster==i]))/sum(cluster==i))/
                  length(m.genes))
        tmp[is.na(tmp)] = 0
      }
      final[i,] = as.numeric(tmp)
    }
    final = final/rowSums(final,na.rm = TRUE)
    final.2 = final
    m = apply(final.2,1,function(x) (max(x)-min(x))/2)
    for (i in z){
      final.2[i,][final.2[i,]<m[i]] = 0
      c.type = names(which(final.2[i,]!=0))
      if (length(c.type)==1){
        if (verbose==TRUE){
          cat(paste(c.names[i], ": ", c.type, sep=""),fill=TRUE)
        }
        class[[i]] = c.type
      }
      if (length(c.type)==2){
        if (verbose==TRUE){
          cat(paste(c.names[i], ": ", c.type[1],"/",c.type[2],sep=""),fill=TRUE)
        }
        class[[i]] = paste(c.type[1],"/",c.type[2],sep="")
      }
      if (length(c.type)>2){
        if (verbose==TRUE){
          cat(paste(c.names[i], ": Unconclusive (more than 2 cell types
                    attributed to the cluster)",sep=""),fill=TRUE)
        }
        class[[i]] = "Undefined"
      }
    }
    names(class) = c.names
    pheatmap(t(final.2),cluster_cols = FALSE)
    types = do.call(rbind,class)
    if (write==TRUE){
      fwrite(data.frame(final),paste("./cluster-analysis/cluster_types.txt",
                                     sep=""),sep="\t")
    }
  }
  res = list(diff.genes,class,types,final.2)
  names(res) = c("diff.genes","class","types","matrix")
  return(res)
}
