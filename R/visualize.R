#' @title Visualize interactions
#' @description Creates chord diagrams from the interactions tables.
#'
#' @details `show.in` gives the elements of `signal` to be displayed in the plot
#' window.
#' @details
#' `write.in` gives the elements of `signal` to be written as pdf files
#' in the *images* folder.
#' @details
#' If `write.out` is TRUE, then the function writes a
#' pdf file with a summary of the all the interactions of `signal` as a chord
#' diagram.
#' @details
#' `limit` is the maximum number of interactions displayed on one chord
#' diagram. Raising this limit over 30 may decrease the visibility.
#'
#' @param signal a list of data frames result of the **cell_signaling()**
#' function
#' @param show.in a vector of which elements of ```signal``` must be shown
#' @param write.in a vector of which elements of ```signal``` must be written
#' @param write.out a logical
#' @param method a string (usually relative to the experiment)
#' @param limit a value between 1 and number of interactions
#'
#' @return The function returns images in the plot window of Rstudio and images
#' in the pdf format in the *images* folder.
#' @export
#' @importFrom circlize chordDiagramFromDataFrame
#' @importFrom circlize circos.trackPlotRegion
#' @importFrom circlize get.cell.meta.data
#' @importFrom circlize circos.text
#' @importFrom circlize circos.axis
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices heat.colors
#' @importFrom grDevices pdf
#' @importFrom graphics mtext
#' @importFrom graphics par
#' @importFrom graphics points
#' @importFrom graphics polygon
#' @importFrom graphics segments
#' @importFrom graphics text
#'
#' @examples
#' int.1 <- matrix(c("gene 1","gene 1", "gene 2", "gene 3"),ncol=2)
#' colnames(int.1) <- c("cluster 1","cluster 2" )
#' int.2 <- matrix(c("gene 1","gene 4","gene 4","gene 2","gene 3","gene 3"),
#' ncol=2)
#' colnames(int.2) <- c("cluster 1","cluster 3" )
#' signal <- list(int.1,int.2)
#' names(signal) <- c("1-2","1-3")
#' visualize_interactions(signal)
visualize_interactions <- function(signal,show.in=NULL,write.in=NULL,write.out=FALSE,
                     method="default",limit=30){
  options(warn=-1)
  if (dir.exists("images")==FALSE & (is.null(write.in)==FALSE |
                                     write.out==TRUE)){
    dir.create("images")
  }
  c.names <- NULL
  for (i in seq_len(length(signal))){
    c.names=c(c.names,colnames(signal[[i]])[seq_len(2)])
  }
  c.names <- unique(c.names)
  cols <- c(rainbow(length(c.names)))
  names(cols) <- c.names
  opar <- par()

  ## Chord diagram of interactions between cluster types -------------

  nn <- NULL
  s <- NULL
  for (i in seq_len(length(signal))){
    if (is.null(dim(signal[[i]]))==TRUE){
      nn <- rbind(nn,names(signal[[i]])[seq_len(2)])
    } else {
      nn <- rbind(nn,colnames(signal[[i]])[seq_len(2)])
    }
    s <- c(s,nrow(signal[[i]]))
  }
  tmp <- cbind(nn,as.numeric(s))
  tmp <- tmp[order(as.numeric(tmp[,3])),]
  name <- tmp[,1]
  feature <- tmp[,2]
  score <- sort(as.numeric(tmp[,3]))
  dat <- data.frame(name,feature)
  cr <- colorRampPalette(c("#FF9900","#CC0000"))(max(score)-min(score)+1)
  if (write.out==TRUE){
    pdf(paste("./images/ChordDiagram",method,"all.pdf",sep="_"))
    chordDiagramFromDataFrame(dat, annotationTrack="grid",
                              preAllocateTracks=1, directional=1,
                              direction.type="arrows",
                              link.arr.length=0.15,
                              link.arr.width=0.15,
                              link.arr.type="triangle",
                              link.arr.lty=par("lty"),
                              link.arr.lwd=par("lwd"),
                              link.arr.col="black",
                              grid.col=cols,col=cr[score-min(score)+1],
                              big.gap=5, small.gap=2)
    circos.trackPlotRegion(track.index=1, panel.fun=function(x, y) {
      xlim=get.cell.meta.data("xlim")
      ylim=get.cell.meta.data("ylim")
      sector.name=get.cell.meta.data("sector.index")
      circos.text(mean(xlim), ylim[1] + .1, sector.name, facing="inside",
                  niceFacing=FALSE, adj=c(0.5, 0),cex=1.5)
      circos.axis(h="top",labels=FALSE, major.tick=FALSE,
                  major.tick.length=1, sector.index=sector.name,
                  track.index=2)
    }, bg.border=NA)
    bx <- par("usr")
    coords <- c(bx[1]*0.95,bx[1]*0.82,bx[4]*0.94,bx[4]*0.53)
    l <- length(cr)
    dy <- (coords[3]-coords[4])/l
    for (i in seq_len(l)){
      x <- c(coords[1],coords[2],coords[2],coords[1])
      y <- c(coords[4]+dy*(i-1),coords[4]+dy*(i-1),coords[4]+dy*i,coords[4]+dy*i)
      polygon(x,y,col=cr[i],border=cr[i])
    }
    text(coords[1]*0.925,coords[3]*0.97,labels=as.character(max(score)),
         col="white")
    text(coords[1]*0.92,coords[4]*1.05,labels=as.character(min(score)),
         col="white")
    text((coords[1]+coords[2])/2,coords[3]*1.03,labels="number")
    dev.off()
  }
  chordDiagramFromDataFrame(dat, annotationTrack="grid",
                            preAllocateTracks=1,
                            directional=1,
                            direction.type="arrows",
                            link.arr.length=0.15,
                            link.arr.width=0.15,
                            link.arr.type="triangle",
                            link.arr.lty=par("lty"),
                            link.arr.lwd=par("lwd"), link.arr.col="black",
                            grid.col=cols,col=cr[score-min(score)+1],
                            big.gap=5, small.gap=2)
  circos.trackPlotRegion(track.index=1, panel.fun=function(x, y) {
    xlim=get.cell.meta.data("xlim")
    ylim=get.cell.meta.data("ylim")
    sector.name=get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, facing="inside",
                niceFacing=FALSE, adj=c(0.5, 0),cex=1.5)
    circos.axis(h="top",labels=FALSE, major.tick=FALSE,
                major.tick.length=1, sector.index=sector.name,
                track.index=2)
  }, bg.border=NA)
  bx <- par("usr")
  coords <- c(bx[1]*0.95,bx[1]*0.82,bx[4]*0.94,bx[4]*0.53)
  l <- length(cr)
  dy <- (coords[3]-coords[4])/l
  for (i in seq_len(l)){
    x <- c(coords[1],coords[2],coords[2],coords[1])
    y <- c(coords[4]+dy*(i-1),coords[4]+dy*(i-1),coords[4]+dy*i,coords[4]+dy*i)
    polygon(x,y,col=cr[i],border=cr[i])
  }
  text(coords[1]*0.925,coords[3]*0.97,labels=as.character(max(score)),
       col="white")
  text(coords[1]*0.92,coords[4]*1.05,labels=as.character(min(score)),
       col="white")
  text((coords[1]+coords[2])/2,coords[3]*1.03,labels="number")


  ## Chord diagram of lig/rec interactions --------------

  for (n in write.in){
    if (is.null(dim(signal[[n]]))==TRUE){
      cnames <- names(signal[[n]])
      tmp <- matrix(signal[[n]],ncol=3)
    } else {
      cnames <- colnames(signal[[n]])
      tmp <- signal[[n]]
      tmp <- tmp[order(as.numeric(tmp[,4]),decreasing=TRUE),]
    }
    if (nrow(tmp)>limit){
      tmp <- tmp[seq_len(limit),]
    }
    tmp <- tmp[order(as.numeric(tmp[,4])),]
    name <- tmp[,2]
    feature <- tmp[,1]
    score <- as.numeric(tmp[,4])*100
    cr <- colorRampPalette(c("lightblue","magenta"))(max(score)-min(score)+1)
    score <- score - min(score) + 1
    dat <- data.frame(name,feature)
    cluster.col <- as.character(c(rep(cols[cnames[2]],length(unique(tmp[,2]))),
                                 rep(cols[cnames[1]],length(unique(tmp[,1])))))
    names(cluster.col) <- c(unique(tmp[,2]),unique(tmp[,1]))
    link.col <- rep("dodgerblue3",nrow(dat))
    link.col[tmp[,3]=="specific"] <- "gray25"
    link.lwd <- rep(1,nrow(dat))
    link.lwd[tmp[,3]=="specific"] <- 3
    link.width <- rep(0.12,nrow(dat))
    link.width[tmp[,3]=="specific"] <- 0.15
    pdf(paste("./images/",cnames[1],"-",cnames[2],"_",
              method,"_diagram.pdf",sep=""))
    chordDiagramFromDataFrame(dat, annotationTrack="grid",
                              preAllocateTracks=1,
                              directional=-1,
                              direction.type="arrows",
                              link.arr.length=link.width,
                              link.arr.width=link.width,
                              link.arr.type="circle",
                              link.arr.lty=par("lty"),
                              link.arr.lwd=link.lwd, link.arr.col=link.col,
                              grid.col=cluster.col,col=cr[score],big.gap=5,
                              small.gap=0.7)
    circos.trackPlotRegion(track.index=1, panel.fun=function(x, y) {
      xlim=get.cell.meta.data("xlim")
      ylim=get.cell.meta.data("ylim")
      sector.name=get.cell.meta.data("sector.index")
      circos.text(mean(xlim), ylim[1] + .1, sector.name, facing="clockwise",
                  niceFacing=TRUE, adj=c(0, 0.5))
      circos.axis(h="top",labels=FALSE,minor.ticks=FALSE,
                  major.tick.length=2,
                  major.at=c(xlim), sector.index=sector.name,
                  track.index=2)
    }, bg.border=NA)
    bx <- par("usr")
    coords <- c(bx[1]*0.95,bx[1]*0.82,bx[4]*0.94,bx[4]*0.53)
    l <- length(cr)
    dy <- (coords[3]-coords[4])/l
    for (i in seq_len(l)){
      x <- c(coords[1],coords[2],coords[2],coords[1])
      y <- c(coords[4]+dy*(i-1),coords[4]+dy*(i-1),coords[4]+dy*i,coords[4]+dy*i)
      polygon(x,y,col=cr[i],border=cr[i])
    }
    text(coords[1]*0.925,coords[3]*0.97,labels=as.character(round(max(
      as.numeric(tmp[,4]))*100)/100),col="white")
    text(coords[1]*0.925,coords[4]*1.05,labels=as.character(round(min(
      as.numeric(tmp[,4]))*100)/100),col="white")
    text((coords[1]+coords[2])/2,coords[3]*1.03,labels="score")
    legend("topright",fill=rev(unique(cluster.col)),legend=cnames[seq_len(2)])
    dev.off()
  }
  for (n in show.in){
    if (is.null(dim(signal[[n]]))==TRUE){
      cnames <- names(signal[[n]])
      tmp <- matrix(signal[[n]],ncol=3)
    } else {
      cnames <- colnames(signal[[n]])
      tmp <- signal[[n]]
      tmp <- tmp[order(as.numeric(tmp[,4]),decreasing=TRUE),]
    }
    if (nrow(tmp)>limit){
      tmp <- tmp[seq_len(limit),]
    }
    tmp <- tmp[order(as.numeric(tmp[,4]),decreasing=FALSE),]
    name <- tmp[,2]
    feature <- tmp[,1]
    score <- as.numeric(tmp[,4])*100
    cr <- colorRampPalette(c("lightblue","magenta"))(max(score)-min(score)+1)
    score <- score - min(score) + 1
    dat <- data.frame(name,feature)
    cluster.col <- as.character(c(rep(cols[cnames[2]],length(unique(tmp[,2]))),
                                 rep(cols[cnames[1]],length(unique(tmp[,1])))))
    names(cluster.col) <- c(unique(tmp[,2]),unique(tmp[,1]))
    link.col <- rep("dodgerblue3",nrow(dat))
    link.col[tmp[,3]=="specific"] <- "gray25"
    link.lwd <- rep(1,nrow(dat))
    link.lwd[tmp[,3]=="specific"] <- 3
    link.width <- rep(0.12,nrow(dat))
    link.width[tmp[,3]=="specific"] <- 0.15
    chordDiagramFromDataFrame(dat, annotationTrack="grid",
                              preAllocateTracks=1,
                              directional=-1,
                              direction.type="arrows",
                              link.arr.length=link.width,
                              link.arr.width=link.width,
                              link.arr.type="circle",
                              link.arr.lty=par("lty"),
                              link.arr.lwd=link.lwd, link.arr.col=link.col,
                              grid.col=cluster.col,col=cr[score],big.gap=5,
                              small.gap=0.2)
    circos.trackPlotRegion(track.index=1, panel.fun=function(x, y) {
      xlim=get.cell.meta.data("xlim")
      ylim=get.cell.meta.data("ylim")
      sector.name=get.cell.meta.data("sector.index")
      circos.text(mean(xlim), ylim[1] + .1, sector.name, facing="clockwise",
                  niceFacing=TRUE, adj=c(0, 0.5),
                  cex=ifelse(nrow(dat)>50, yes=0.5,no=1))
      circos.axis(h="top",labels=FALSE,minor.ticks=FALSE,
                  major.at=c(xlim), major.tick.length=2,
                  sector.index=sector.name, track.index=2)

    }, bg.border=NA)
    bx <- par("usr")
    coords <- c(bx[1]*0.95,bx[1]*0.82,bx[4]*0.94,bx[4]*0.53)
    l <- length(cr)
    dy <- (coords[3]-coords[4])/l
    for (i in seq_len(l)){
      x <- c(coords[1],coords[2],coords[2],coords[1])
      y <- c(coords[4]+dy*(i-1),coords[4]+dy*(i-1),coords[4]+dy*i,coords[4]+dy*i)
      polygon(x,y,col=cr[i],border=cr[i])
    }
    text(coords[1]*0.925,coords[3]*0.97,
         labels=as.character(round(max(as.numeric(tmp[,4]))*100)/100),
         col="white")
    text(coords[1]*0.925,coords[4]*1.05,
         labels=as.character(round(min(as.numeric(tmp[,4]))*100)/100),
         col="white")
    text((coords[1]+coords[2])/2,coords[3]*1.03,labels="score")
    legend("topright",fill=rev(unique(cluster.col)),legend=cnames[seq_len(2)])
  }

  par(opar)
}
