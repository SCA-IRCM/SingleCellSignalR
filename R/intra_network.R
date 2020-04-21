#' @title intra network
#' @description Computes intracellular networks linked to genes of interest.
#'
#' @details `signal` is a list containing the cell-cell interaction tables. It
#' is the result of the **cell_signaling()** function.
#' @details
#' `cell.prop` is set to 0.2
#' by default to avoid unreadable downstream networks. However if the calculated
#' network is too small or non-existent (or too big) the user can try lower
#' (or higher) values.
#' @details
#' If the user does not set `c.names`, the clusters will be
#' named from 1 to the maximum number of clusters (cluster 1, cluster 2, ...).
#' The user can exploit the `c.names` vector in the list returned by the
#' **cell_classifier()** function for this purpose. The user can also provide
#' her own cluster names.
#' @details
#' `species` must be equal to "homo sapiens" or
#' "mus musculus". In the case of mouse data, the function converts mouse genes
#' in human orthologs (according to Ensembl) such that the Reactome/KEGG
#' interaction database can be exploited, and finally output genes are converted
#' back to mouse.
#' @details
#' If `write` is TRUE, then the function writes two different
#' files. A graphML file in the *network* folder for intracellular
#' interactions downstream the gene of interest (goi) named
#' "intracell_network_coi-receptors.graphml". A text file in the *network*
#' folder containing the information about the pathways in which the interactions
#' are in, named "intracell_network_pathway_analysis_coi-receptors.txt".
#'
#' @param goi gene of interest (typically a receptor)
#' @param data a data frame of n rows (genes) and m columns (cells) of read or
#' UMI counts (note : rownames(data)=genes)
#' @param genes a character vector of HUGO official gene symbols of length n
#' @param cluster a numeric vector of length m
#' @param coi name of the cluster of interest
#' @param cell.prop a threshold, only the genes expressed in this proportion of
#' the cells of the coi will be taken into account
#' @param c.names (optional) cluster names
#' @param signal (optional) a list (result of the **cell_signaling()** function)
#' @param species "homo sapiens" or "mus musculus"
#' @param write a logical (if TRUE writes graphML and text files for the
#' internal networks)
#' @param plot a logical
#' @param add.lig a logical (if TRUE adds the goi associated ligands from
#' signal to the network)
#' @param connected a logical (if TRUE keeps only the genes connected to
#' the goi)
#' @param verbose a logical
#'
#' @return The function returns a list containing the internal networks
#' linked to the genes of interest (goi)
#'
#' @export
#'
#' @importFrom foreach foreach %do%
#' @importFrom stats phyper
#' @importFrom multtest mt.rawp2adjp
#' @importFrom utils write.table
#' @import igraph
#'
#' @examples
#' data <- matrix(runif(1000,0,1),nrow=5,ncol=200)
#' genes <- c("A2M","LRP1","AANAT","MTNR1A","ACE")
#' cluster <- c(rep(1,100),rep(2,100))
#' intra_network(goi=c("TGFBR1","ERBB2"),data,genes,cluster,coi="T-cells")
intra_network <- function(goi,data,genes,cluster,coi,cell.prop=0.2,c.names=NULL,
                         signal=NULL,write=TRUE,plot=TRUE,add.lig=TRUE,
                         species=c("homo sapiens","mus musculus"),
                         connected=FALSE,verbose=TRUE){
  if (dir.exists("networks")==FALSE & write==TRUE){
    dir.create("networks")
  }
  if (is.null(c.names)==TRUE){
    c.names <- paste("cluster",seq_len(max(cluster)))
  }
  if (min(cluster)!=1){
    cluster <- cluster + 1 - min(cluster)
  }
  if (length(c.names)!=max(cluster) | sum(duplicated(c.names))>0 |
      grepl("/",paste(c.names,collapse =""))){
    stop("The length of c.names must be equal to the number of clusters and must
        contain no duplicates. The cluster names must not include special
        characters")
  }
  if (!is.element(coi,c.names)){
    stop(paste(coi,"must be included in c.names.","If c.names is not provided, it is set to cluster 1, cluster 2, ...,
        cluster N. WIth N the maximum number of clusters"))
  }
  opar <- par()
  species <- match.arg(species)
  if (species=='mus musculus'){
    Hs2mm <- mm2Hs[,1]
    mm2Hs <- mm2Hs[,2]
    names(mm2Hs) <- Hs2mm
    names(Hs2mm) <- as.character(mm2Hs)
    m.names <- mm2Hs[rownames(data)]
    data <- subset(data,(!is.na(m.names)))
    m.names <- m.names[!is.na(m.names)]
    rownames(data) <- as.character(m.names)
    goi <- mm2Hs[goi]
  }
  pw.names <- strsplit(PwC_ReactomeKEGG$pathway,";")
  pw.sizes <- table(unlist(pw.names))
  max.pw.size <- 500
  good.pw <- pw.sizes[pw.sizes<=max.pw.size]

  data.tmp <- data[,cluster==which(c.names==coi)]
  data.tmp <- data.tmp[rowSums(data.tmp)>0,]

  # expressed genes
  good <- apply(data.tmp,1,function(x) sum(x>0)/ncol(data.tmp)>cell.prop)
  visible.genes <- unique(c(rownames(data.tmp)[good],goi))
  visible.n <- PwC_ReactomeKEGG[PwC_ReactomeKEGG$a.gn%in%visible.genes &
                                 PwC_ReactomeKEGG$b.gn%in%visible.genes,]
  red.visible.n <- simplify_interactions(visible.n,LRdb)
  res <- list()
  qq <- 0

  for (receptors in goi){
    qq <- qq+1
    if (!is.element(receptors,rownames(data.tmp))){
      cat(paste0(receptors," is not expressed in ",coi),fill=TRUE)
    } else {
      # receptor containing pathways ------------
      contains.receptors <- intersect(unlist(
        pw.names[PwC_ReactomeKEGG$a.gn%in%receptors |
                   PwC_ReactomeKEGG$b.gn%in%receptors]),names(good.pw))
      if (verbose == TRUE){
        cat(paste0("Patwhay(s) that include ",receptors,":"),fill=TRUE)
        for (i in contains.receptors){
          cat(paste0("   - ",i),fill=TRUE)
        }
      }
      if (sum(grepl("added",contains.receptors))>1){
        contains.receptors <- contains.receptors[contains.receptors!="added"]
      }
      contain.n=NULL
      for (i in contains.receptors){
        contain.n <- rbind(contain.n,PwC_ReactomeKEGG[
          grepl(i,PwC_ReactomeKEGG$pathway),])
      }
      contain.n <- unique(rbind(contain.n,PwC_ReactomeKEGG[
        PwC_ReactomeKEGG$a.gn%in%visible.genes &
          PwC_ReactomeKEGG$b.gn%in%receptors,],
        PwC_ReactomeKEGG[PwC_ReactomeKEGG$a.gn%in%receptors &
                           PwC_ReactomeKEGG$b.gn%in%visible.genes,]))
      red.contain.n <- simplify_interactions(contain.n)

      # intersect corr network and receptor-containing network ------------
      key.visible <- paste(red.visible.n$a.gn,red.visible.n$b.gn,sep="|")
      key.contain <- paste(red.contain.n$a.gn,red.contain.n$b.gn,sep="|")
      net.n <- red.visible.n[key.visible%in%key.contain,]

      if (nrow(net.n)>0){
        # add ligands
        add.net <- NULL
        nam=NULL
        siz=NULL
        if (is.null(signal)==FALSE & add.lig==TRUE){
          for (i in names(signal)){
            if (grepl(paste0(coi,"$"),i)){
              tmp=signal[[i]]
              if (species=="mus musculus"){
                m.1 <- mm2Hs[tmp[,1]]
                m.2 <- mm2Hs[tmp[,2]]
                tmp <- subset(tmp,(!is.na(m.1)))
                tmp <- subset(tmp,(!is.na(m.2)))
                m.1 <- m.1[!is.na(m.1)]
                m.2 <- m.2[!is.na(m.2)]
                tmp[,1] <- m.1
                tmp[,2] <- m.2
              }
              if (is.element(receptors,tmp[,2])){
                siz <- c(siz,rowMeans(data[tmp[tmp[,2]%in%receptors,1],
                                          cluster==as.numeric(which(
                                            c.names%in%colnames(tmp)[1]))]))
                nam.tmp <- rep(colnames(tmp)[1],nrow(tmp[tmp[,2]%in%receptors,]))
                colnames(tmp)[seq_len(2)] <- c("a.gn","b.gn")
                tmp[,1] <- paste(nam.tmp[1],tmp[,1],sep="-")
                add.net <- rbind(add.net,tmp[tmp[,2]%in%receptors,])
                nam <- c(nam,nam.tmp)
              }
            }
          }
          if (is.null(add.net)==FALSE){
            add.net <- cbind(add.net[,seq_len(2)],nam,location="extra",
                            type="control")
          }
        }
        net.tmp <- cbind(net.n[,seq_len(2)],nam=rep(coi,nrow(net.n)),
                        location="intra",type=net.n$type)
        net.f <- rbind(add.net,net.tmp)
        if (species=="mus musculus"){
          net.f[,1] <- Hs2mm[net.f[,1]]
          net.f[,2] <- Hs2mm[net.f[,2]]
        }
        g.net <- graph_from_data_frame(net.f,directed=TRUE)
        g.net.tmp <- as.undirected(g.net)

        y <- shortest_paths(g.net.tmp,receptors,V(g.net.tmp))
        g.net.tmp <- graph_from_data_frame(net.f[,seq_len(2)],directed=FALSE)
        V(g.net.tmp)$status <- "pw.related"
        V(g.net.tmp)$status[
          which(unique(c(net.f$a.gn,net.f$b.gn)) %in% receptors)] =
          "gene.of.interest"
        V(g.net.tmp)$status[which(unique(c(net.f$a.gn,net.f$b.gn)) %in%
                                    net.f$a.gn[net.f$location=="extra"])] =
          "ligand"
        E(g.net.tmp)$int.type <- as.character(net.f$type)

        y <- unlist(lapply(y$vpath,function(x) length(x)))
        y <- y-1
        if (sum(y==-1)>0 & connected==FALSE){
          y[y==-1] <- 1
        }
        if (sum(y==-1)>0 & connected==TRUE){
          nam.tmp <- unique(c(net.f$a.gn,net.f$b.gn))[y==-1]
          net.f <- net.f[!net.f$a.gn%in%nam.tmp & !net.f$b.gn%in%nam.tmp,]
          y <- y[y!=-1]
          g.net <- graph_from_data_frame(net.f,directed=TRUE)
        }
        names(y) <- unique(c(net.f$a.gn,net.f$b.gn))
        y[net.f$a.gn[net.f$location=="extra"]]=-1
        x=vector("numeric",length=length(y))
        names(x) <- unique(c(net.f$a.gn,net.f$b.gn))
        if (sum(net.f$location=="extra")){
          x[net.f$a.gn[net.f$location=="extra"]] <- (seq(0,length(
            net.f$a.gn[net.f$location=="extra"])*4-1,4)-max(seq(0,length(
              net.f$a.gn[net.f$location=="extra"])*4-1,4))/2)
        }
        x[receptors] <- 0
        for (j in seq_len(max(y))){
          if (sum(y==j)>1){
            x[y==j] <- seq(-10,10,19/sum(y==j))[seq_len(sum(y==j))]
          }
        }
        l <- cbind(x,-y)
        for (i in seq_len(max(y))){
          if(sum(y==i)>1){
            index <- which(y==i)
            l[index[seq(1,length(index),2)],2] <-
              l[index[seq(1,length(index),2)],2]-0.2
            l[index[seq(2,length(index),2)],2] <-
              l[index[seq(2,length(index),2)],2]+0.2
          }
        }

        rownames(l) <- unique(c(net.f$a.gn,net.f$b.gn))

        if (sum(net.f$location %in% "intra")==0){
          V(g.net)$size[y==-1] <- (log(siz)+5)*4
          V(g.net)$size[V(g.net)$size<exp(-5)] <- exp(-5)
        } else {
          V(g.net)$size <- (log(c(siz,rowMeans(data.tmp[
            unique(c(net.n$a.gn,net.n$b.gn)),])))+5)*4
          V(g.net)$size[V(g.net)$size<exp(-5)]=exp(-5)
        }
        V(g.net)$size[which(unique(c(net.f$a.gn,net.f$b.gn)) %in% receptors)] <-
          35
        V(g.net)$shape <- c("circle")
        V(g.net)$shape[which(unique(c(net.f$a.gn,net.f$b.gn)) %in% receptors)] <-
          c("rectangle")
        V(g.net)$color <- c("lightcyan")
        V(g.net)$color[which(unique(c(net.f$a.gn,net.f$b.gn)) %in% receptors)] <-
          c("indianred1")
        V(g.net)$color[which(unique(c(net.f$a.gn,net.f$b.gn)) %in%
                               net.f$a.gn[net.f$location=="extra"])] <-
          c("lightyellow1")
        V(g.net)$label.color <- "black"
        V(g.net)$label.dist[which(!unique(c(net.f$a.gn,net.f$b.gn))
                                  %in% receptors)] <- 1
        for (i in unique(l[,2])){
          if (i!=0){
            V(g.net)$label.degree[l[,2]==i] <- rep(c(pi/2,-pi/2),length(y))[
              seq_len(sum(l[,2]==i))]
          }
        }
        V(g.net)$label.dist[which(unique(c(net.f$a.gn,net.f$b.gn)) %in%
                                    receptors)] <- 0
        V(g.net)$label.degree[which(unique(c(net.f$a.gn,net.f$b.gn)) %in%
                                      receptors)] <- 0
        E(g.net)$arrow.mode <- as.numeric(grepl("control",net.f$type))*2
        E(g.net)$arrow.size <- rep(0.4,nrow(net.f))
        E(g.net)$color <- "gray"

        # pathway enrichment statistics
        g.pw.names <- strsplit(net.n$pathway,";")
        g.pw.names <- g.pw.names[g.pw.names!="added"]
        g.pw.sizes <- table(unlist(g.pw.names))
        N <- nrow(PwC_ReactomeKEGG)
        n <- nrow(net.n)
        pw.table <- foreach(pw=names(g.pw.sizes),.combine=rbind) %do% {
          K <- pw.sizes[pw]
          if (K <= max.pw.size){
            k <- g.pw.sizes[pw]
            pval <- 1-phyper(q=k-1,m=K,n=N-K,k=n)
            data.frame(pathway=pw,in.pw=k,pw.size=K,pval=pval,
                       stringsAsFactors=FALSE)
          }
          else
            NULL
        }

        par(mfrow=c(2,1))
        par(las=2)
        par(mar=c(0,0,1,0))
        plot(g.net, layout=l,main=coi)

        if (is.null(pw.table)==FALSE){
          rawp <- pw.table$pval
          if (length(rawp)==1){
            adj <- rawp
            qval <- adj
          } else {
            adj <- mt.rawp2adjp(rawp,"BH")
            qval <- adj$adjp[order(adj$index),"BH"]
          }

          pw.table <- cbind(pw.table,data.frame(qval=qval))
          rownames(pw.table) <- NULL
          pw.table <- pw.table[pw.table$qval<0.05,]
          if (nrow(pw.table)>0){
            pw.table <- pw.table[order(pw.table$in.pw,decreasing=FALSE),]
            nc <- max(nchar(as.character(pw.table$pathway)))*0.3

            mtext(paste(receptors,"related pathways"),side=1,line=2,
                  las=FALSE)
            par(mar=c(2,nc,4,2))
            barplot((pw.table$in.pw),horiz=TRUE,
                    names.arg=(paste(pw.table$pathway,"*")),
                    cex.names=0.7,col="azure2",border="gray60")
          } else {
            if (verbose==TRUE){
              cat("No significant associated pathway",fill=TRUE)
            }
          }
        } else {
          if (verbose ==TRUE){
            cat(paste("No associated genes downstream",receptors, "in",coi),
                fill=TRUE)
          }
        }

        if (write==TRUE){
          if (is.null(pw.table)==FALSE){
            write.table(pw.table[order(pw.table$qval),],
                        file=paste0(
                          './networks/intracell_network_pathway_analysis_',
                                    coi,'-',receptors,'.txt'),sep="\t",
                        quote=FALSE,row.names=FALSE)
          }
          write.graph(g.net.tmp,paste0('./networks/intracell_network_',coi,
                                       '-',receptors,'.graphml'),
                      format="graphml")
        }
        res[[qq]] <- net.f
      }
      else {
        cat("No interactions found.", fill=TRUE)
      }
    }
  }
  par(opar)
  return(res)
}
