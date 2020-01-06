#' @title Cell Signaling
#' @description Computes "autocrine" or "paracrine" interactions between cell
#' clusters.
#'
#' @param data a data frame of n rows (genes) and m columns (cells) of read or
#' UMI counts (note : rownames(data)=genes)
#' @param genes a character vector of HUGO official gene symbols of length n
#' @param cluster a numeric vector of length m
#' @param int.type "autocrine" or "paracrine"
#' @param c.names (optional) cluster names
#' @param s.score LRscore threshold
#' @param logFC a number, the log fold-change threshold for differentially
#' expressed genes
#' @param species "homo sapiens" or "mus musculus"
#' @param tol a tolerance parameter for balancing between "autocrine" and
#' "paracrine" interactions
#' @param write a logical
#' @param verbose a logical
#'
#' @return The function returns "autocrine" or "paracrine" interaction lists.
#'
#' @export
#' @importFrom stats median
#' @importFrom stats p.adjust
#' @importFrom stats median pnorm
#' @importFrom stats median sd
#' @importFrom igraph graph_from_data_frame write.graph
#'
#' @examples
#' data=matrix(runif(1000,0,1),nrow=5,ncol=200)
#' rownames(data) = c("A2M","LRP1","AANAT","MTNR1A","ACE")
#' cluster=c(rep(1,100),rep(2,100))
#' cell_signaling(data,rownames(data),cluster,int.type="paracrine",write=FALSE)
cell_signaling = function(data, genes,
                          cluster,int.type=c("paracrine","autocrine"),
                          c.names=NULL,s.score=0.5,logFC=log2(1.5),
                          species=c("homo sapiens","mus musculus"),
                          tol=0,write=TRUE,verbose=TRUE){
  if (dir.exists("cell-signaling")==FALSE & write==TRUE){
    dir.create("cell-signaling")
  }
  if (is.null(c.names)==TRUE){
    c.names = paste("cluster",seq_len(max(cluster)))
  }
  if (min(cluster)!=1){
    cluster = cluster + 1 - min(cluster)
  }
  if (length(c.names)!=max(cluster) | sum(duplicated(c.names))>0 |
      grepl("/",paste(c.names,collapse =""))){
    cat("The length of c.names must be equal to the number of clusters
        and must contain no duplicates. The cluster names must not include
        special characters",fill=TRUE)
    return()
  }
  int.type = match.arg(int.type)
  species = match.arg(species)

  rownames(data) = genes
  z = seq_len(max(cluster))
  lig = unique(LRdb$ligand)
  rec = unique(LRdb$receptor)
  data = data.frame(data)
  data = data[rowSums(data)>0,]
  med = sum(data)/(nrow(data)*ncol(data))

  if (species=='mus musculus'){
    mm2Hs = mm2Hs[!is.na(mm2Hs$`Mouse orthology confidence`) &
                    mm2Hs$`Mouse orthology confidence`==1,c(3,5)]
    mm2Hs = subset(mm2Hs,!duplicated(mm2Hs$`Gene name`))
    mm2Hs = subset(mm2Hs,mm2Hs$`Mouse gene name`!="a")
    Hs2mm = mm2Hs[,1]
    mm2Hs = mm2Hs[,2]
    names(mm2Hs) = as.character(Hs2mm)
    names(Hs2mm) = as.character(mm2Hs)
    m.names = mm2Hs[rownames(data)]
    data = subset(data,(!is.na(m.names)))
    m.names = m.names[!is.na(m.names)]
    rownames(data)=as.character(m.names)
  }


  ## Autocrine -------------------
  if (int.type=="autocrine"){
    if (verbose==TRUE){
      cat("Autocrine signaling: ",fill=TRUE)
    }
    auto = list()
    k=0
    int=NULL
    n.int=NULL
    if (verbose==TRUE){
      cat("Checking for cell/cell signaling:",fill=TRUE)
    }
    for (i in z){
      tmp = data[,cluster==i]
      tmp = tmp[rowSums(tmp)>0,]
      if (sum(is.element(lig, rownames(tmp)))>0){
        lig.tmp = rownames(tmp)[is.element(rownames(tmp),lig)]
        #lig.tmp = lig.tmp[!is.element(lig.tmp,gene.list[[i]])]
      } else {lig.tmp=NULL}

      final.tmp = LRdb[is.element(LRdb$ligand,lig.tmp),seq_len(2)]
      final.tmp = data.frame(final.tmp,as.character(
        rep("autocrine|paracrine",sum(is.element(LRdb$ligand,lig.tmp)))))
      m.lig = rowSums(tmp[unique(final.tmp[,1]),])/sum(cluster==i)
      names(m.lig) = unique(final.tmp[,1])

      if (sum(is.element(rec, rownames(tmp)))>0){
        rec.tmp = rownames(tmp)[is.element(rownames(tmp[apply(
          tmp,1,function(x) sum(x>0))>tol*ncol(tmp),]),rec)]
      } else {rec.tmp=NULL}

      for (j in z){
        temp = data[,cluster==j]
        temp = temp[rowSums(temp)>0,]
        if (sum(is.element(rec, rownames(temp)))>0){
          rec.temp = rownames(temp)[is.element(rownames(temp),rec)]
          #rec.temp = rec.temp[!is.element(rec.temp,gene.list[[j]])]
        } else {rec.temp=NULL}
        rec.temp = rec.temp[is.element(rec.temp,rec.tmp)]
        m.rec = rowSums(data.frame(temp[rec.temp,]))/sum(cluster==j)
        names(m.rec) = rec.temp

        final = final.tmp[is.element(final.tmp$receptor,rec.temp),]
        final = cbind(final,score(m.lig[final$ligand],m.rec[final$receptor],
                                  med))

        colnames(final) = c(c.names[i],c.names[j],"interaction type","LRscore")

        if (i==j){
          final$`interaction type`="autocrine"
        }

        final = final[final[,4]>s.score,]
        final = final[order(final[,4],decreasing = TRUE),]

        if (species=="mus musculus"){
          final[,1] = Hs2mm[as.character(final[,1])]
          final[,2] = Hs2mm[as.character(final[,2])]
        }

        if (nrow(final)>0){
          k=k+1
          auto[[k]] = final
          if (verbose==TRUE){
            cat(paste(nrow(final),"interactions from",c.names[i],
                      "to",c.names[j]),fill=TRUE)
          }
          int = c(int,paste(i,"-",j,sep=""))
          n.int = c(n.int,paste(c.names[i],"-",c.names[j],sep=""))
          gr = graph_from_data_frame(final,directed=FALSE)
          if (write==TRUE){
            fwrite(data.frame(final),paste("./cell-signaling/LR_interactions_",
                                           c.names[i],"-",c.names[j],"-",
                                           int.type,".txt",sep=""),sep="\t")
          }
        }
      }
    }
    if (k!=0){
      names(auto) = n.int
    }
  }


  ## Paracrine -------------------
  if (int.type=="paracrine"){
    gene.list = vector("list",max(cluster))
    for (i in z){
      if (file.exists(paste("./cluster-analysis/table_dge_",c.names[i],
                            ".txt",sep=""))==TRUE){
        resu = fread(paste("./cluster-analysis/table_dge_",c.names[i],
                           ".txt",sep=""),data.table = FALSE)
        gene.list[[i]] = resu$genes[resu$logFC>logFC]
        if (species == "mus musculus"){
          gene.list[[i]] = mm2Hs[gene.list[[i]]]
        }
        gene.list[[i]] = gene.list[[i]][!is.na(gene.list[[i]])]
      } else {
        gene.list[[i]] = "none"
        cat(paste("No such file as table_dge_",c.names[i],
                  ".txt in the cluster-analysis folder", sep=""),fill=TRUE)
      }
    }

    if (verbose==TRUE){
      cat("Paracrine signaling: ",fill=TRUE)
    }
    para = list()
    k=0
    int=NULL
    n.int=NULL
    if (verbose==TRUE){
      cat("Checking for signaling between cell types",fill=TRUE)
    }
    for (i in z){
      if (sum(cluster==i)>1){
        tmp = data[,cluster==i]
        tmp = tmp[rowSums(tmp)>0,]
        if (sum(is.element(lig, rownames(tmp)))>0){
          lig.tmp = rownames(tmp)[is.element(rownames(tmp),lig)]
          #lig.tmp = lig.tmp[!is.element(lig.tmp,gene.list[[i]])]
        } else {lig.tmp=NULL}

        final.tmp = LRdb[is.element(LRdb$ligand,lig.tmp),seq_len(2)]
        final.tmp = data.frame(final.tmp,as.character(
          rep("paracrine",sum(is.element(LRdb$ligand,lig.tmp)))),
          stringsAsFactors = FALSE)
        m.lig = rowSums(tmp[unique(final.tmp[,1]),])/sum(cluster==i)
        names(m.lig) = unique(final.tmp[,1])

        if (sum(is.element(rec, rownames(tmp)))>0){
          rec.tmp = rownames(tmp)[is.element(rownames(tmp[
            apply(tmp,1,function(x) sum(x>0))>tol*ncol(tmp),]),rec)]
        } else {rec.tmp=NULL}

        for (j in z[-i]){
          if (sum(cluster==j)>1){
            temp = data[,cluster==j]
            temp = temp[rowSums(temp)>0,]
            if (sum(is.element(rec, rownames(temp)))>0){
              rec.temp = rownames(temp)[is.element(rownames(temp),rec)]
              #rec.temp = rec.temp[!is.element(rec.temp,gene.list[[j]])]
            } else {rec.temp=NULL}
            rec.temp = rec.temp[!is.element(rec.temp,rec.tmp)]
            m.rec = rowSums(data.frame(temp[rec.temp,]))/sum(cluster==j)
            names(m.rec) = rec.temp

            final = final.tmp[is.element(final.tmp$receptor,rec.temp),]
            final = cbind(final,score(m.lig[final$ligand],m.rec[final$receptor],
                                      med))
            exclus = final$ligand %in% gene.list[[i]] & final$receptor %in%
              gene.list[[j]]
            if (sum(exclus)!=0){
              f.exclu = final[exclus,]
              final = final[!(final$ligand %in% gene.list[[i]] & final$receptor
                              %in% gene.list[[j]]),]
              f.exclu[,3] = "specific"
              final = rbind(f.exclu,final)
            }

            colnames(final) = c(c.names[i],c.names[j],"interaction type",
                                "LRscore")
            final = final[final[,4]>s.score,]
            final = final[order(final[,4],decreasing = TRUE),]

            if (species=="mus musculus"){
              final[,1] = Hs2mm[as.character(final[,1])]
              final[,2] = Hs2mm[as.character(final[,2])]
            }

            if (nrow(final)>0){
              k=k+1
              para[[k]] = final
              if (verbose==TRUE){
                cat(paste(nrow(final),"interactions from",c.names[i],"to",
                          c.names[j]),fill=TRUE)
              }
              int = c(int,paste(i,"-",j,sep=""))
              n.int = c(n.int,paste(c.names[i],"-",c.names[j],sep=""))
              gr = graph_from_data_frame(final,directed=FALSE)
              if (write==TRUE){
                fwrite(data.frame(final),paste(
                  "./cell-signaling/LR_interactions_",c.names[i],"-",c.names[j],
                  "-",int.type,".txt",sep=""),sep="\t")
              }
            }
          }
        }
      }
    }
    if (k!=0){
      names(para) = n.int
    }
  }


  ## Returns ---------------------
  if (int.type=="autocrine"){
    return(auto)
  }
  if (int.type=="paracrine"){
    return(para)
  }
}

#' Calculation of the LRscore
#'
#' @param l a value (or a vector) of the mean ligand expression in the
#' secreting cluster
#' @param r a value (or a vector) of the mean receptor expression in the
#' receiving cluster
#' @param s a value for scaling the score (usually the mean of the whole
#' read count table, the median or another similar value is possible),
#' must be over 0
#'
#' @return a value or a vector
#' @export
#'
#' @examples
#' l=1
#' r=9
#' s=5
#' score(l,r,s)
score = function(l,r,s){
  L=l^(1/2)
  R=r^(1/2)
  S=s
  sc=L*R/(S+L*R)
  return(sc)
}
