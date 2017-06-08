

dotplot_pipesteamed<-function
(steamed, # result of running a pipe_steam like pipe_dbscan
 chopMethod="tsne",# tsne or PCA
 steamMethod="dbscan",
 cName="group"
 ){

  sampTab<-steamed[['steamed']][['sampTab']]

  if(chopMethod=='tsne'){
    choppedDat<-steamed[['cp_tsne']][['choppedDat']]
  }
  else{
    choppedDat<-steamed[['cp_pca']][['choppedDat']]
  }
  datTab<-choppedDat[rownames(sampTab),]
  colnames(datTab)[1:2]<-c("dim.1", "dim.2")

  xres<-cbind(sampTab, datTab)
  xi<-which(colnames(xres)==cName)
  colnames(xres)[xi]<-"groupX"

  ColorRamp <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(unique(xres$groupX)))
  ggplot(xres, aes(x=dim.1, y=dim.2, colour=groupX) ) + geom_point(pch=19, alpha=3/4, size=1) + theme_bw() + scale_colour_manual(values=ColorRamp) #+ facet_wrap( ~ k, nrow=3)
}



#' heatmap of the classification result
#'
#' Heatmap of the classification result.
#' @param classMat classMat
#' @param isBig is this a big heatmap? TRUE or FALSE
#' @param cluster_cols cluster_cols
#'
#' @return nothing
#'
#' @examples
#' cn_HmClass(cnRes, isBig=TRUE)
#'
#' @export
hmClass<-function
(classMat, 
 isBig=FALSE,
 cluster_cols=FALSE
){
 
  cools<-colorRampPalette(c("black", "limegreen", "yellow"))( 100 )
  bcol<-'white';
  if(isBig){
    bcol<-NA;
  }
  pheatmap(classMat,
    col=cools,
    breaks=seq(from=0, to=1, length.out=100),
    border_color=bcol,
    cluster_rows = FALSE,
    cluster_cols = cluster_cols)
  # classification heatmap
}


#' @export
hmcellsgenes<-function
(expDat,
 geneList,
 sampTab,# must of a colname 'group' corresponding to names of geneList
 n_genes=5,
 limits=c(-3,3),
 fontsize_row=NULL){


   dLevel<-"group"
   
   geneList<-geneList[order(as.numeric(names(geneList)))]
   groupNames<-names(geneList)

   stX<-data.frame()
   for(groupName in groupNames){
    stX<-rbind(stX, sampTab[sampTab$group==groupName,])
   }
   stX<-droplevels(stX)

   topx<-unlist(lapply(geneList, length))
   topx[which(topx>=n_genes)]<-n_genes
#   topx<-topx[order(as.numeric(names(topx)))]
   geneList<-geneList[ names(topx) ]
   genes<-list()
   for(groupName in groupNames){
     if(topx[groupName]>0){
       genes[[groupName]]<-geneList[[groupName]][1:as.vector(topx[groupName])]
     }
   }
   genes<-unlist(genes)
   cat("genes: ", length(genes),"\n") 
   stX<-stX[order(stX[,dLevel]),]
  
   expDat<-expDat[,rownames(stX)]
   value <- t(scale(t(expDat[genes,])))
   value[value < limits[1]] <- limits[1]
   value[value > limits[2]] <- limits[2]


   gene_annotation<-data.frame(group=rep(names(topx), topx))
   genes2<-make.unique(genes)
   rownames(gene_annotation)<-genes2
   rownames(value)<-genes2
  
   cell_annotation<-data.frame(group=stX[,dLevel])
   rownames(cell_annotation)<-rownames(stX)
 

   

    xcol <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(groupNames))
    names(xcol) <- groupNames
    anno_colors <- list(group = xcol)

    pheatmap(value, cluster_rows = FALSE, cluster_cols = FALSE,
        show_colnames = F, annotation_names_row = FALSE,
        annotation_col = cell_annotation, annotation_row = gene_annotation,
        annotation_names_col = FALSE, annotation_colors = anno_colors, fontsize_row=fontsize_row);
   list(c1=cell_annotation, g1=gene_annotation)
}


#' @export
hmgenesSimple<-function
(expDat,
  cRow=FALSE,
  cCol=FALSE,
 limits=c(-3,3),
 toScale=FALSE,
 fontsize_row=NULL,
 anc=FALSE){
  

  value<-expDat
  if(toScale){
    value <- t(scale(t(expDat)))
  }

  
  value[value < limits[1]] <- limits[1]
  value[value > limits[2]] <- limits[2]

    pheatmap(value, 
      cluster_rows = cRow, 
      cluster_cols = cCol,
        show_colnames = anc, annotation_names_row = FALSE, scale='none',
#        clustering_distance_rows='euclidean',
        clustering_distance_rows='correlation',
        clustering_distance_cols="correlation", 
       #clustering_distance_cols="euclidean", 
        clustering_method="average",
        annotation_names_col = FALSE, fontsize_row=fontsize_row)
  #  value;
}

#' @export
hm_varGenes<-function
(washed,
 steamed, # from steam_pipe so has cp_pca
  cName="group",
  cRow=FALSE,
  cCol=FALSE,
  limits=c(-3,3),
  toScale=FALSE,
  fontsize_row=NULL){
  
  expDat<-washed[['expDat']]
  varGenes<-steamed[['cp_pca']][['varGenes']]
  sampTab<-steamed[['steamed']][['sampTab']]
  expDat<-expDat[varGenes,rownames(sampTab)]
  
  annTab<-data.frame(group=sampTab[,cName])
  rownames(annTab)<-rownames(sampTab)

  value<-expDat
  if(toScale){
    value <- t(scale(t(expDat)))
  }
  
  value[value < limits[1]] <- limits[1]
  value[value > limits[2]] <- limits[2]

 
  groupNames<-unique(annTab[,cName])
  xcol <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(groupNames))
  names(xcol) <- groupNames
  anno_colors <- list(group = xcol)

  xx<-data.frame(group=as.factor(annTab[,cName]))
  rownames(xx)<-rownames(annTab)

  pheatmap(value, cluster_rows = cRow, cluster_cols = cCol,
        show_colnames = FALSE, annotation_names_row = FALSE,
##        annotation_col = annTab,
 ##clustering_distance_rows='correlation',
        annotation_col = xx,
        annotation_names_col = FALSE, annotation_colors = anno_colors, fontsize_row=fontsize_row)
}



#' @export
hmgenesSimple2<-function
(expDat,
  annTab,
  cName="group",
  cRow=FALSE,
  cCol=FALSE,
 limits=c(-3,3),
 toScale=FALSE,
 fontsize_row=NULL){
  

  value<-expDat
  if(toScale){
    value <- t(scale(t(expDat)))
  }

  
  value[value < limits[1]] <- limits[1]
  value[value > limits[2]] <- limits[2]

 
  groupNames<-unique(annTab[,cName])
  xcol <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(groupNames))
  names(xcol) <- groupNames
  anno_colors <- list(group = xcol)

if(FALSE){
  pheatmap(value, 
    cluster_rows = cRow, 
    cluster_cols = cCol,
    show_colnames = F, annotation_names_row = FALSE, annotation_names_col = TRUE, scale='none',
#        clustering_distance_rows='euclidean',
        clustering_distance_rows='correlation',
       clustering_distance_cols="correlation", 
       #clustering_distance_cols="euclidean", 
        clustering_method="average",
        annotation_col = annTab, annotation_colors = anno_colors, fontsize_row=fontsize_row)
  #  value;
}

  xx<-data.frame(group=as.factor(annTab[,cName]))
  rownames(xx)<-rownames(annTab)

  pheatmap(value, cluster_rows = cRow, cluster_cols = cCol,
        show_colnames = FALSE, annotation_names_row = FALSE,
##        annotation_col = annTab,
        annotation_col = xx,
        annotation_names_col = FALSE, annotation_colors = anno_colors, fontsize_row=fontsize_row)
}


#' plot tsne results
#'
#' plot tsne results
#'
#' @param mcRes, which must have TSNE.1 and TSNE.2, and group columns
#' @param legend whether to display it
#'
#' @return ggplot
#' 
#' @export
#'
plotDBscan<-function(mcRes, legend=FALSE)
{
  ColorRamp <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(unique(mcRes$group)))
  if(legend){
    ans<-ggplot(mcRes, aes(x=TSNE.1, y=TSNE.2, colour=group) ) + geom_point(pch=19, alpha=3/4, size=.5) + theme_bw() + scale_colour_manual(values=ColorRamp)
  }
  else{
    ans<-ggplot(mcRes, aes(x=TSNE.1, y=TSNE.2, colour=group) ) + geom_point(pch=19, alpha=3/4, size=.5) + theme_bw() + scale_colour_manual(values=ColorRamp) + theme(legend.position="none")
  }
  ans
}

#' @export
sc_plot_statTab<-function# multi plot of mu, alpha, mean, cov, fano, and mean vs cov
(statTab){


  plot1<-ggplot(statTab, aes(x=mu)) + geom_histogram(aes(y=..density..), colour="black", fill="#238b45") + theme_bw() + ggtitle("mu")
  plot2<-ggplot(statTab, aes(x=alpha)) + geom_histogram(aes(y=..density..), colour="black", fill="#88419d") + theme_bw() + ggtitle("alpha")
  plot3<-ggplot(statTab, aes(x=overall_mean)) + geom_histogram(aes(y=..density..), colour="black", fill="#2b8cbe") + theme_bw() + ggtitle("overall_mean")
  plot4<-ggplot(statTab, aes(x=fano)) + geom_histogram(aes(y=..density..), colour="black", fill="#d7301f") + theme_bw() + ggtitle("fano")
  plot5<-ggplot(statTab, aes(x=cov)) + geom_histogram(aes(y=..density..), colour="black", fill="#feb24c") + theme_bw() + ggtitle("cov")
  plot6<-ggplot(statTab, aes(x=overall_mean, y=sd)) + geom_point(pch=19, alpha=2/4, size=1, colour="#9e9ac8") + theme_bw() + ggtitle("overall mean vs sd")

  multiplot(plot1, plot2, plot3, plot4, plot5, plot6, cols=2)

}

#' @export
ptsne<-function(xres, cname="study_id")
{
  xi<-which(colnames(xres)==cname)
  colnames(xres)[xi]<-"group"
  ColorRamp <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(unique(xres$group)))
  ggplot(xres, aes(x=TSNE.1, y=TSNE.2, colour=group) ) + geom_point(pch=19, alpha=3/4, size=1) + theme_bw() + scale_colour_manual(values=ColorRamp) #+ facet_wrap( ~ k, nrow=3)
}

#' @export
tsneClass<-function
(classRes, #result of butter_classify()
 steamed #from a pipe_wash like pipe_dbscan()
 ){
  sampTab<-steamed[['steamed']][['sampTab']]
  choppedDat<-steamed[['cp_tsne']][['choppedDat']]
  datTab<-choppedDat[rownames(sampTab),]
  datTab<-cbind(datTab, t(classRes[,rownames(sampTab)]))
  classNames<-rownames(classRes)
  datTab<-as.data.frame(datTab)
  tsneMult(datTab, classNames)
}

#' @export
tsneMult<-function# facet tsne plot by gene
(tsneDat, # cols TSNE.1, TNSE.2, and genes
 genesToPlot, # genes to plot
 colorPal="BuPu",
 revCol=TRUE
 ){
  require(tidyr)
  tsne<-tsneDat[,c("TSNE.1", "TSNE.2", genesToPlot)]
  tsneLong<-gather_(tsne, "gene", "expression", genesToPlot)

  if(revCol){
    ColorRamp <- rev(colorRampPalette(rev(brewer.pal(n = 7,name = colorPal)))(100))[10:100]
  }
  else{
    ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7,name = colorPal)))(100)
  }
  ggplot(tsneLong, aes(x=TSNE.1, y=TSNE.2, colour=expression) ) + 
  geom_point(pch=19, alpha=2/4, size=.25) + 
  theme_bw() + 
  scale_colour_gradientn(colours=ColorRamp) +
  facet_wrap( ~ gene)
 
}


# generic red/blue heatmap
#' @export
mp_hmVars<-function# basic heatmap
(expDat,
 ### numeric matrix
 genes,
 ### rownames to include in the heatmap
 main='',
 ### optional title for teh top of the HM
 clusterR=T,
 ### whether to cluster the Rows
 clusterC=F,
 ### whether to cluster the columns
 scale='row',
 ### normalize the data: 'row', 'column', or 'none'
 big=FALSE,
 ### if not big then add cell separators
 dist=utils_myDist,
 ###
 margin=c(12,6),
 ###
 RowSideColors=NULL,
 ###
 ColSideColors=NULL,
 ###
 cexCol=1,
 cexRow=1,
 ccol=''
){
  
  genes<-intersect(rownames(expDat), genes);
  
  if(length(ccol)==4){
    ccol<-colorpanel(ccol$n, ccol$low,ccol$mid, ccol$high);
  }
  else{
    ccol<-bluered(100);
  }
  if(is.null(RowSideColors)){
    RowSideColors<-rep('white', length(genes));
  }
  if(is.null(ColSideColors)){
    ColSideColors<-rep('white', ncol(expDat));
  }
  if(!clusterR & ! clusterC){
    dendrogram='none';
  }
  if(clusterR & clusterC){
    dendrogram='both';
  }
  if(clusterR & !clusterC){
    dendrogram='row';
  }
  if(!clusterR & clusterC){
    dendrogram='column';
  }
  if(!big){
    heatmap.2(expDat[genes,],
              #col=bluered(100),
              col=ccol,
              scale=scale,
              trace='none',
              key=T,
              dendrogram=dendrogram,
              Rowv=clusterR,
              Colv=clusterC,
              density.info='none',
              margin=margin,
              colsep=seq(ncol(expDat)),
              rowsep=seq(length(genes)),
              sepcol='white',
              sepwidth=c(0.001,0.00),
              main=main,
              dist=dist,
              RowSideColors=RowSideColors,
              ColSideColors=ColSideColors,
              cexCol=cexCol,
              cexRow=cexRow);
  }
  else{
    heatmap.2(expDat[genes,],col=ccol, scale=scale, trace='none', key=T,dendrogram=dendrogram,Rowv=clusterR,Colv=clusterC,density.info='none',margin=margin,main=main,dist=dist,labRow='',labCol='',RowSideColors=RowSideColors,ColSideColors=ColSideColors,cexCol=cexCol,cexRow=cexRow);
  }
}




# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
# See: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#' @export
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



