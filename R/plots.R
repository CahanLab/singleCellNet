#'
reorderCells<-function
(grpList){
  curr_grps<-grpList[[1]]
  curr_cells<-names(curr_grps)

  for(i in 2:length(grpList)){
    prior_grps<-curr_grps
    prior_cells<-curr_cells

    xList<-grpList[[i]]
    prior_names<-sort(unique(prior_grps))
    for(j in seq(length(prior_names))){
      pname<-prior_names[j]
      xi<-which(prior_grps==pname)
      tmpCells<-prior_cells[xi]
      tmpRes<-sort(xList[tmpCells])

      curr_cells[xi]<-names(tmpRes)
      curr_grps[xi] <-tmpRes
    }
  }
  names(curr_grps)<-curr_cells
  curr_grps
}






#' @export
corplot_sub<-function
(gpaRes,
 expDat,
 prop=0.1,
 pSide=FALSE,
 minCount=20,
 llevel=1){

  ###orderedCells<-reorderCells(gpaRes$grp_list)
  orderedCells<-reorderCells(gpaRes$grp_list)

  ssamp<-subsample_min(orderedCells, prop=prop, minCount=minCount)
  ##ssamp<-subsample_min(gpaRes$groups, prop=prop, minCount=minCount)
###  ssamp<-subsample_min(gpaRes$grp_list[[2]], prop=prop, minCount=minCount)



  genes<-getVarFromList(gpaRes, llevel=llevel)
  
###
###  xcor<-cor(expDat[genes,names(ssamp)])
###

  xcor<- dist(t(expDat[genes,names(ssamp)]))

  ##xcor<-as.matrix(gpaRes$results[["L1_G1"]]$gpRes$xdist)[names(ssamp), names(ssamp)]
  maxX<-max(xcor)
  xcor<-as.matrix((maxX - xcor)/maxX)

  llevels<-length(gpaRes$grp_list)

  xx<-gpaRes$grp_list[[2]][names(ssamp)]
  xx<-data.frame(level_1=as.factor(xx))
  cnames<-c("level_1")

  if(llevels>2){
    for(i in 3:llevels){
      cnames<-append(cnames, paste0("level_", i-1))
      danss<-gpaRes$grp_list[[i]][names(ssamp)]
      xx<-cbind(xx, as.factor(danss))
    }
    colnames(xx)<-cnames
  }
 
  if(pSide){
    #topGenes<-gpaRes$groupTree$Get("topGenes")


    topGenes<-getTopGenesList(gpaRes$results[["L1_G1"]],7)


    ### xy<-data.frame(levelX=xx[,ncol(xx)], genes=rep("", nrow(xx)))
    xy<-data.frame(levelX=xx[,1], genes=rep("", nrow(xx)))

    rLabels<-rep("", nrow(xy))
    grpNames<-unique(xy$levelX)
    for(grpName in grpNames){
      cat("***",grpName,"***\n")
      xi<-which(xy$levelX==grpName)
      coord<-ceiling( (max(xi)-min(xi)) / 2 ) + min(xi)
      cat(grpName,"  xi:",xi[1], "length: ", length(xi), "coord: ", coord,"\n")
      rLabels[coord]<-topGenes[grpName]
    } 


   ## ColorRamp <- colorRampPalette(brewer.pal(n = 9,name = "Blues"))(100)

    rownames(xy)<-rownames(xx)
    pheatmap(xcor,
     ## color=ColorRamp,
      cluster_rows = FALSE, 
      cluster_cols = FALSE,  
      show_colnames = FALSE,
      show_rownames=TRUE,
      annotation_names_row = FALSE, 
      annotation_col = xx,
#      annotation_row = xy,
      labels_row=rLabels,
      fontsize_row=5)
  }
  else{
    pheatmap(xcor,
      cluster_rows = FALSE, 
      cluster_cols = FALSE,  
      show_colnames = FALSE,
      show_rownames=FALSE,
      annotation_names_row = FALSE, 
      annotation_col = xx)
    }

}


#' @export
corplot_triangle<-function
(gpaRes,
 expDat,
 prop=0.1,
 pSide=FALSE,
 minCount=20,
 llevel=1){

  ###orderedCells<-reorderCells(gpaRes$grp_list)
  orderedCells<-reorderCells(gpaRes$grp_list)

  ssamp<-subsample_min(orderedCells, prop=prop, minCount=minCount)
  ##ssamp<-subsample_min(gpaRes$groups, prop=prop, minCount=minCount)
###  ssamp<-subsample_min(gpaRes$grp_list[[2]], prop=prop, minCount=minCount)



  genes<-getVarFromList(gpaRes, llevel=llevel)
  
###
###  xcor<-cor(expDat[genes,names(ssamp)])
###

###  xcor<- dist(t(expDat[genes,names(ssamp)]))

  xcor<-as.matrix(gpaRes$results[["L1_G1"]]$gpRes$xdist)[names(ssamp), names(ssamp)]
  maxX<-max(xcor)
  xcor<-as.matrix((maxX - xcor)/maxX)

 xcor[lower.tri(xcor)]<-0

  llevels<-length(gpaRes$grp_list)

  xx<-gpaRes$grp_list[[2]][names(ssamp)]
  xx<-data.frame(level_1=as.factor(xx))
  cnames<-c("level_1")

  if(llevels>2){
    for(i in 3:llevels){
      cnames<-append(cnames, paste0("level_", i-1))
      danss<-gpaRes$grp_list[[i]][names(ssamp)]
      xx<-cbind(xx, as.factor(danss))
    }
    colnames(xx)<-cnames
  }
 
  if(pSide){
    #topGenes<-gpaRes$groupTree$Get("topGenes")


    topGenes<-getTopGenesList(gpaRes$results[["L1_G1"]],7)


    ### xy<-data.frame(levelX=xx[,ncol(xx)], genes=rep("", nrow(xx)))
    xy<-data.frame(levelX=xx[,1], genes=rep("", nrow(xx)))

    rLabels<-rep("", nrow(xy))
    grpNames<-unique(xy$levelX)
    for(grpName in grpNames){
      cat("***",grpName,"***\n")
      xi<-which(xy$levelX==grpName)
      coord<-ceiling( (max(xi)-min(xi)) / 2 ) + min(xi)
      cat(grpName,"  xi:",xi[1], "length: ", length(xi), "coord: ", coord,"\n")
      rLabels[coord]<-topGenes[grpName]
    } 


    ColorRamp <- colorRampPalette(brewer.pal(n = 9,name = "Blues"))(100)

    rownames(xy)<-rownames(xx)
    pheatmap(xcor,
      color=ColorRamp,
      cluster_rows = FALSE, 
      cluster_cols = FALSE,  
      show_colnames = FALSE,
      show_rownames=TRUE,
      annotation_names_row = FALSE, 
      annotation_col = xx,
#      annotation_row = xy,
      labels_row=rLabels,
      fontsize_row=5)
  }
  else{
    pheatmap(xcor,
      cluster_rows = FALSE, 
      cluster_cols = FALSE,  
      show_colnames = FALSE,
      show_rownames=FALSE,
      annotation_names_row = FALSE, 
      annotation_col = xx)
    }

}


getVarFromList<-function(
  gpaRes,
  llevel=1
  ){

   gresNames<-names(gpaRes$results)
   charMatch<-paste0("L",llevel,"_")

   gresNames<-gresNames[grep(charMatch, gresNames)]


if(FALSE){
   if(llevel==1){
    tmpAns<-gpaRes$results[[1]]$gpRes$pcaRes$varGenes
   }
   else{
      tmpAns<-vector()
      for(i in seq(length(gpaRes$results))){
       tmpAns<-append(tmpAns, gpaRes$results[[i]]$gpRes$pcaRes$varGenes)
      }
    } 
  }
  tmpAns<-vector()
  for(gName in gresNames){
    tmpAns<-append(tmpAns, gpaRes$results[[gName]]$gpRes$pcaRes$varGenes)
  }
  sort(unique(tmpAns))
}


subsample_min<-function
(namedVect,
  prop=0.10,
  minCount=20){
  groups<-unique(namedVect)
  grpCounts<-table(namedVect)
  newVect<-vector()
  for(grp in groups){
    cat(grp, ": ", ceiling(prop*grpCounts[grp]),"\n")
    numToSamp<-max(minCount, ceiling(prop*grpCounts[grp]))
    xi<-which(namedVect==grp)
    tmpVect<-sample(namedVect[xi], numToSamp)
    newVect<-append(newVect, tmpVect)
  }
  newVect
}




#' @export
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
hm_genes<-function
(washed,
 steamed, # from steam_pipe so has cp_pca
 genes,
  cName="group",
  cRow=FALSE,
  cCol=FALSE,
  limits=c(-3,3),
  toScale=FALSE,
  fontsize_row=NULL){
  
  expDat<-washed[['expDat']]
  ###varGenes<-steamed[['cp_pca']][['varGenes']]
  sampTab<-steamed[['steamed']][['sampTab']]
  expDat<-expDat[genes,rownames(sampTab)]
  
  annTab<-data.frame(group=sampTab[,cName])
  rownames(annTab)<-rownames(sampTab)

  value<-expDat
  if(toScale){
    value <- t(scale(t(expDat)))
  }
  
  value[value < limits[1]] <- limits[1]
  value[value > limits[2]] <- limits[2]

 
 ### groupNames<-unique(annTab[,cName])
  groupNames<-unique(sampTab[,cName])

  xcol <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(groupNames))
  names(xcol) <- groupNames
  anno_colors <- list(group = xcol)

  xx<-data.frame(group=as.factor(annTab[,"group"]))
  rownames(xx)<-rownames(annTab)

  pheatmap(value, cluster_rows = cRow, cluster_cols = cCol,
        show_colnames = FALSE, annotation_names_row = FALSE,
##        annotation_col = annTab,
   clustering_distance_rows='correlation',
   clustering_distance_cols='correlation',
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


#' plot gpa res
#'
#' plot gpa res
#'
#' @param gpaRes,
#' @param legend whether to display it
#'
#' @return ggplot
#' 
#' @export
#'
plotGPALevel<-function(gpaRes, llevel="L1_G1",legend=FALSE)
{

  xdat<-gpaRes$results[[llevel]]$gpRes$pcaRes$pcaRes$x
  grps<-gpaRes$results[[llevel]]$bundleRes$result
  xdat<-xdat[names(grps),]
  aDat<-data.frame(pc1=xdat[,1], pc2=xdat[,2], group=grps) 
 
  ColorRamp <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(unique(aDat$group)))
  if(legend){
    ans<-ggplot(aDat, aes(x=pc1, y=pc2, colour=group) ) + geom_point(pch=19, alpha=3/4, size=.5) + theme_bw() + scale_colour_manual(values=ColorRamp)
  }
  else{
    ans<-ggplot(aDat, aes(x=pc1, y=pc2, colour=group) ) + geom_point(pch=19, alpha=3/4, size=.5) + theme_bw() + scale_colour_manual(values=ColorRamp) + theme(legend.position="none")
  }
  ans
}


#' plot gpa res
#'
#' plot gpa res
#'
#' @param gpaRes,
#' @param legend whether to display it
#'
#' @return ggplot
#' 
#' @export
#'
plotGPA1<-function(gpaRes, legend=FALSE)
{

  aDat<-data.frame(pc1=gpaRes$pcaRes$pcaRes$x[,1], pc2=gpaRes$pcaRes$pcaRes$x[,2], group=gpaRes$groups) 
  ColorRamp <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(unique(aDat$group)))
  if(legend){
    ans<-ggplot(aDat, aes(x=pc1, y=pc2, colour=group) ) + geom_point(pch=19, alpha=3/4, size=.5) + theme_bw() + scale_colour_manual(values=ColorRamp)
  }
  else{
    ans<-ggplot(aDat, aes(x=pc1, y=pc2, colour=group) ) + geom_point(pch=19, alpha=3/4, size=.5) + theme_bw() + scale_colour_manual(values=ColorRamp) + theme(legend.position="none")
  }
  ans
}

#' plot gpa_recurse res
#'
#' plot gpa_recurse res
#'
#' @param gpaRes,
#' @param legend whether to display it
#'
#' @return ggplot
#' 
#' @export
#'
plotGPArec<-function(gpaRes, pcaLevel=1, grpLevel=1,legend=FALSE)
{
 
  if(grpLevel>1){
    prev<-gpaRes$grp_list[[grpLevel-1]]
  }
  else{
    prev<-rep("1", length(gpaRes$grp_list[[grpLevel]]))
  }
  aDat<-data.frame(pc1=gpaRes$results[[pcaLevel]]$pcs[,1], pc2=gpaRes$results[[pcaLevel]]$pcs[,2], group=gpaRes$grp_list[[grpLevel]]) 
  ColorRamp <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(unique(aDat$group)))
  if(legend){
    ans<-ggplot(aDat, aes(x=pc1, y=pc2, colour=group, shape=prev) ) + geom_point( alpha=3/4, size=.75) + theme_bw() + scale_colour_manual(values=ColorRamp)
  }
  else{
    ans<-ggplot(aDat, aes(x=pc1, y=pc2, colour=group, shape=prev) ) + geom_point(pch=19, alpha=3/4, size=.75) + theme_bw() + scale_colour_manual(values=ColorRamp) + theme(legend.position="none")
  }
  ans
}

#' plot gpa_recurse res transform coordinates
#'
#' plot gpa_recurse res
#'
#' @param gpaRes,
#' @param legend whether to display it
#'
#' @return ggplot
#' 
#' @export
#'
plotGPArecTrans<-function(gpaRes, pcaLevel=1, grpLevel=1,legend=FALSE)
{
 
  pc_l1<-gpaRes$results[[1]]$pcs[,1:2]

  runningPCs<-pc_l1
  ci<-2

  while(ci<=pcaLevel){
    weight<-(pcaLevel+1-ci)/pcaLevel
    cat("CI ",ci, " weight: ",weight,"\n")
    tmpPCs<-weight*gpaRes$results[[ci]]$pcs[,1:2]
    runningPCs<-tmpPCs+runningPCs
    ci<-ci+1
  }

  if(grpLevel>1){
    prev<-gpaRes$grp_list[[grpLevel-1]]
  }
  else{
    prev<-rep("1", length(gpaRes$grp_list[[grpLevel]]))
  }
  aDat<-data.frame(pc1=runningPCs[,1], pc2=runningPCs[,2], group=gpaRes$grp_list[[grpLevel]]) 
  ColorRamp <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(unique(aDat$group)))
  if(legend){
    ans<-ggplot(aDat, aes(x=pc1, y=pc2, colour=group, shape=prev) ) + geom_point( alpha=3/4, size=.75) + theme_bw() + scale_colour_manual(values=ColorRamp)
  }
  else{
    ans<-ggplot(aDat, aes(x=pc1, y=pc2, colour=group, shape=prev) ) + geom_point( alpha=3/4, size=.75) + theme_bw() + scale_colour_manual(values=ColorRamp) + theme(legend.position="none")
  }
  ans

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

simple_tnse<-function
(tsneRes){
  xres<-as.data.frame(tsneRes)
  ggplot(xres, aes(x=TSNE.1, y=TSNE.2) ) + geom_point(pch=19, alpha=2/4, size=1) + theme_bw() 
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

tsneGenes<-function
(washed,
 steamed,
 genes)
{
  sampTab<-steamed[['steamed']][['sampTab']]
  choppedDat<-steamed[['cp_tsne']][['choppedDat']]
  expDat<-washed$expDat[genes,rownames(sampTab)]

  datTab<-choppedDat[rownames(sampTab),]
  datTab<-cbind(datTab, t(expDat))
  datTab<-as.data.frame(datTab)
  tsneMult(datTab, genes)
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

#' @export
tsneMultsimp<-function# facet tsne plot by gene
(tsneDat, # cols TSNE.1, TNSE.2
  expDat,
 genesToPlot, # genes to plot
 colorPal="BuPu",
 revCol=TRUE,
 toScale=TRUE,
 limits=c(-3,3)
 ){

  require(tidyr)
  value<-expDat[genesToPlot,]

  if(toScale){
    value <- t(scale(t(value)))
  }
  value[value < limits[1]] <- limits[1]
  value[value > limits[2]] <- limits[2]


  tsne<-as.data.frame(tsneDat[,c("TSNE.1", "TSNE.2")])
  tsne<-cbind(tsne, t(value))

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

#' plot gpa res
#'
#' plot gpa res
#'
#' @param gpaRes,
#' @param title
#' @param legend whether to display it
#'
#' @return ggplot
#' 
#' @export
#'
plotGPA<-function(gpaRes, title='', legend=FALSE)
{

  aDat<-data.frame(pc1=gpaRes$pcaRes$pcaRes$x[,1], pc2=gpaRes$pcaRes$pcaRes$x[,2], group=gpaRes$groups) 
  ColorRamp <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(unique(aDat$group)))
  if(legend){
    ans<-ggplot(aDat, aes(x=pc1, y=pc2, colour=group) ) + geom_point(pch=19, alpha=3/4, size=.5) + theme_bw() + scale_colour_manual(values=ColorRamp)
  }
  else{
    ans<-ggplot(aDat, aes(x=pc1, y=pc2, colour=group) ) + geom_point(pch=19, alpha=3/4, size=.5) + theme_bw() + scale_colour_manual(values=ColorRamp) + theme(legend.position="none")
  }
  ans + ggtitle(title)
}

plotGPAann<-function(gpaRes, sampTab, dLevel="prefix", title='')
{

  
  aDat<-data.frame(pc1=gpaRes$pcaRes$pcaRes$x[,1], pc2=gpaRes$pcaRes$pcaRes$x[,2], group=gpaRes$groups) 
  rownames(aDat)<-rownames(gpaRes$pcs)
  aDat<-cbind(aDat, ann=sampTab[rownames(aDat),dLevel])
  ColorRamp <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(unique(aDat$group)))
  
    ans<-ggplot(aDat, aes(x=pc1, y=pc2, colour=group, shape=ann) ) + geom_point(alpha=3/4, size=.5) + theme_bw() + scale_colour_manual(values=ColorRamp)
  
  ans + ggtitle(title)
}



gpa_multiPlot_ann<-function(
 recRes,
 sampTab,
 dLevel="prefix",
 str_level=1
 ){

  pList<-list()
  p1<-plotGPAann(recRes$results[[str_level]]$res[[1]], sampTab=sampTab, dLevel=dLevel, title="Top level")
  xlen<-length(recRes$results[[str_level+1]]$res)
  pList[[1]]<-p1
  for(xi in 1:xlen){
    ttitle<-names(recRes$results[[str_level+1]]$res)[xi]
    cat(ttitle,"\n")
    xtmp<-recRes$results[[str_level+1]]$res[[xi]]
    if(length(xtmp)>1){
      pList[[xi+1]]<-plotGPAann(xtmp, sampTab=sampTab, dLevel=dLevel,title=ttitle)
    }
  }
  multiplot(plotlist=pList, cols=xlen+1)
}


gpa_multiPlot<-function(
 recRes,
 str_level=1
 ){

  pList<-list()
  p1<-plotGPA(recRes$results[[str_level]]$res[[1]], legend=T, title="Top level")
  xlen<-length(recRes$results[[str_level+1]]$res)
  pList[[1]]<-p1
  for(xi in 1:xlen){
    ttitle<-names(recRes$results[[str_level+1]]$res)[xi]
    cat(ttitle,"\n")
    xtmp<-recRes$results[[str_level+1]]$res[[xi]]
    if(length(xtmp)>1){
      pList[[xi+1]]<-plotGPA(xtmp, legend=F, title=ttitle)
    }
  }
  multiplot(plotlist=pList, cols=xlen+1)
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



