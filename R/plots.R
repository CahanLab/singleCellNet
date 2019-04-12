# Patrick Cahan (C) 2017
# patrick.cahan@gmail.com


#' @export
sc_violinClass<-function
(sampTab,
 classRes,
 cellIDCol = "cell_name",
 dLevel="cluster",
 addRand=0,
 threshold=0.20, 
 ncol =1,
 sub_cluster = NA){

  rownames(sampTab) = sampTab[,cellIDCol]
  sids <- rownames(sampTab)
  colnames(sampTab)[which(colnames(sampTab) == dLevel)] = "cluster"
  dLevel = "cluster"
  classRes<-classRes[,sids]
  stQ2<-cbind(sampTab[sids,], t(classRes[,sids]))

  maxX <-apply(classRes, 1, max)
  meaVar <-names(which(maxX>threshold))

  test <- melt(stQ2, id.vars = c(cellIDCol, dLevel), measure.vars =  meaVar)

  cnames <- colnames(test)
  cnames[which(cnames=='value')] <- "classification_score"
  cnames[which(cnames=='variable')] <- "cell_type"
  colnames(test) <- cnames
  xcol = length(unique(sampTab[,dLevel]))
  getPalette <- colorRampPalette(brewer.pal(xcol, "Set2"))
  if(!is.na(sub_cluster)){
    test = test[test$cluster %in% sub_cluster,]
  }

ggplot(test, aes(x = cluster, y = classification_score, fill = cluster)) + ylim(0,1) + geom_violin(scale='width', position='dodge', trim=FALSE) + 
  facet_wrap(~ cell_type, ncol=ncol) + scale_y_continuous(
   # expand = c(0, 0),
     name = "Class score",
     breaks = c(0,  0.50,  1.0),
     labels = c("0", "0.50", "1.0"),
     limits = c(0,1)
   ) +
   coord_cartesian(clip = "off") +
   theme_dviz_hgrid() +
   theme(
     axis.line.x = element_blank(),
     axis.ticks.x = element_blank(), 
     axis.title.y = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    axis.text.x = element_text(size = 8, angle = 45, hjust=1),
    axis.text.y = element_text(size = 8),
    legend.title=element_text(size=10), 
    legend.text=element_text(size=10),
    strip.text.x = element_text(size = 8)
   ) + geom_boxplot(width=0.1, outlier.size=.25) + scale_fill_manual(values = getPalette(xcol))

}


#' @export
prep_umap_class<-function
(classRes,
 sampTab,
 nrand,
 dLevel,
 uCates='all',
 sid='sample_name',
 topPC=10){


  stTmp<-addToST(sampTab, nrand=nrand, sid=sid, dLevels=dLevel)
  stTmp<-assign_cate(classRes, stTmp)
  colnames(stTmp)[2] <- "group"
  if(uCates!='all'){
    uCates<-unique(stTmp[,"category"])
   }
   else{
    ###uCates <- unique(as.vector(sampTab[,dLevel]))
    uCates<-rownames(classRes)
   }
  cat("PCA\n")
  pcRes<-prcomp(t(classRes[uCates,]))
  if(topPC>length(uCates)){
    topPC <- length(uCates)
  }
  cat("UMAP\n")
  uRes<-umap(pcRes$x[,1:topPC], min_dist=.5)
  cat("done")
  stTmp<-cbind(stTmp, uRes$layout)
  colnames(stTmp)[4]<-"umap.1"
  colnames(stTmp)[5]<-"umap.2"
  stTmp
}

#' @export
plot_umap<-function
(preRes){

  ColorRamp <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(unique(preRes$category)))
  ggplot(preRes, aes(x=umap.1, y=umap.2, colour=category) ) + geom_point(pch=19, alpha=3/4, size=1) + theme_bw() + scale_colour_manual(values=ColorRamp) 

}

#' @export
addToST<-function(sampTab, nrand, sid="sample_name", dLevels=c("description1")){

  grpRand<-rep("rand", nrand)
  names(grpRand)<-paste("rand_", 1:nrand, sep='')

  rTab<-data.frame(sid=names(grpRand))
  for(dlev in dLevels){
    rTab<-cbind(rTab, grpRand)
  }
  colnames(rTab)<-c("sid", dLevels)
  qTab<-sampTab[,c(sid, dLevels)]
  colnames(qTab)[1]<-"sid"
  rbind(qTab, rTab)
}



#' @export
assign_cate<-function(classRes, sampTab, cThresh=0){
  topCats<-rownames(classRes)[apply(classRes, 2, which.max)]
  sampTab<-cbind(sampTab, category=topCats)
  sampTab
}

#' @export
plot_attr<-function(classRes, sampTab, nrand, dLevel, sid="sample_name", sub_cluster = NA){
  if(nrand>0){
    stTmp<-addToST(sampTab, nrand=nrand, sid=sid, dLevels=dLevel)} else{
      stTmp = sampTab
    }
  stTmp<-assign_cate(classRes, stTmp)
  colnames(stTmp)[which(colnames(stTmp) == dLevel)]="group"
  getPalette = colorRampPalette(brewer.pal(12, "Paired"))
  myPal = getPalette(length(unique(stTmp$category)))
  p1 = ggplot(stTmp, aes(x=group, fill=category)) +  geom_bar(position = "fill", width=.6) + scale_y_continuous(labels = scales::percent) + scale_fill_manual(values=myPal) +  theme_bw() + coord_flip()
if(is.na(sub_cluster)){
  p1
}else{
  ggplot(stTmp[which(stTmp$group %in% sub_cluster),], aes(x = group, fill = category)) + geom_bar(position = "fill", width = 0.6) + scale_y_continuous(labels = scales::percent) + 
    scale_fill_manual(values = myPal) + theme_bw()+ theme(legend.position = "bottom",,axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
}}


#' @export
makeColList<-function(sampTab,dLevel,cellName="cell_name",palName="Blues"){
    grps<-as.vector(sampTab[,dLevel])
    names(grps)<-as.vector(sampTab[,cellName])
    cells<-names(grps)  
    groupNames<-unique(grps)

    xcol <- colorRampPalette(rev(brewer.pal(n = 9,name = palName)))(length(groupNames)+1)[1:length(groupNames)]
    names(xcol) <- groupNames
    anno_colors <- list(group = xcol)

    xx<-data.frame(group=as.factor(grps))
    rownames(xx)<-cells
    list(ann_col = xx, ann_colors = anno_colors)
}

#' @export
hm_sel_col<-function(
  expDat,
  genes,
  clusCols,# from makeColLis
  maxPerGrp=100,
  cRow=FALSE,
  cCol=FALSE,
  limits=c(0,10),
  toScale=FALSE,
  fontsize_row=4){

  allgenes<-rownames(expDat)
  missingGenes<-setdiff(genes, allgenes)
  if(length(missingGenes)>0){
    cat("Missing genes: ", paste0(missingGenes, collapse=","), "\n")
    genes<-intersect(genes, allgenes)
  }

  value<-expDat[genes,]
  if(toScale){
      value <- t(scale(t(value)))
    }

  value[value < limits[1]] <- limits[1]
  value[value > limits[2]] <- limits[2]

  anno_colors <- clusCols[["ann_colors"]]
  xx <- clusCols[["ann_col"]]
  groupNames<-names(anno_colors$group)
  cells<-colnames(expDat)

  cells2<-vector()
  for(groupName in groupNames){
    xi<-which(grps==groupName)
    if(length(xi)>maxPerGrp){
      tmpCells<-sample(cells[xi], maxPerGrp)
    }
    else{
      tmpCells<-cells[xi]
    }
    cells2<-append(cells2, tmpCells)
  }
  value<-value[,cells2]

   val_col <- colorRampPalette(rev(brewer.pal(n = 12,name = "Spectral")))(25)
  pheatmap(value, cluster_rows = cRow, cluster_cols = cCol, color=val_col,
        show_colnames = FALSE, annotation_names_row = FALSE,
        annotation_col = xx,
        annotation_names_col = FALSE, annotation_colors = anno_colors, fontsize_row=fontsize_row)
}





#' Skyline waterfall
#'
#' Skyline waterfall of classification result
#' @param classRes
#' @param cellType
#' @param sampTab
#' @param dLevel
#' @param yval
#'
#' @return nothing
#'
#' @examples
#' skylineClass(classRes, "blood",sampTab, "timepoint", 0.25)
#'
#' @export
skylineClass<-function(
  classRes, 
  cellType,
  sampTab, 
  dLevel,
  yval,
  cellIdLab="cell_name",
  reorder=TRUE){



  cellVals<-as.vector(sampTab[,dLevel])
  cellAnns<-unique(cellVals)

  # assign cols somehow
  xdf<-cbind(sampTab, t(classRes))

  newXdf<-data.frame();

  # re-order by cellAnn then by classification of cellType
  if(reorder){
    for(group in cellAnns){
      xDat<-xdf[which(xdf[,dLevel]==group),]
     xDat<-xDat[ order(xDat[,cellType]),]
     newXdf<-rbind(newXdf, xDat)
    }
  }
  else{
    newXdf<-xdf
  }
  

  xi<-which(colnames(newXdf) == cellType)
  colnames(newXdf)[xi] <- "vals"

  xi<-which(colnames(newXdf) == dLevel)
  colnames(newXdf)[xi] <- "group"

  newXdf[,cellIdLab]<-factor(newXdf[,cellIdLab], as.vector(newXdf[,cellIdLab]))
  xxi<-which(colnames(newXdf)==cellIdLab)
  colnames(newXdf)[xxi]<-"cell_name"
  #ggplot(data=newXdf, aes(x=cell_name, y=vals, fill=cellAnns)) + geom_bar(stat="identity") + theme_bw()
  ggplot(data=newXdf, aes(x=cell_name, y=vals, fill=group)) + geom_bar(stat="identity", width = 1) + theme_bw() + geom_hline(aes(yintercept=yval), colour="#990000", linetype="dashed") + scale_x_discrete(breaks=NULL) + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + scale_fill_viridis(discrete=TRUE) + ylim(c(0,1))
  #newXdf
}



#' heatmap of the enrichment result
#'
#' Heatmap of tthe enrichment result
#' @param esList returned from ks.extract.more
#' @param threshold threshold 0.05 (holm-corrected pval)
#'
#' @return nothing
#'
#' @examples
#' hm_enr(esList)
#'
#' @export
hm_enr<-function
(esList, # has ES and pval
 threshold=0.05,
 cRows=FALSE, 
 cCols=FALSE,
 fsr=4
){
 
  
  cools <- colorRampPalette(rev(brewer.pal(n = 11,name = "RdBu")))(100)
 # cools<-colorRampPalette(c("blue", "white", "yellow"))( 100 )

  pes <- esList[['ES']]
  pp  <- esList[['pval']]

  pes[which(pp>threshold)]<-0
  mmin<-min(pes)
  mmax<-max(pes)
  mmax<-max(abs(mmin), abs(mmax))


    bcol<-"white"

  pheatmap(pes,
    col=cools,
    breaks=seq(from= -1 * mmax, to=mmax, length.out=100),
    border_color=bcol,
    cluster_rows = cRows,
    cluster_cols = cCols,
    fontsize_row=fsr)
}


#' re-order cells for plotting
#'
#' re-order cells for plotting
#'
#' @param nvect named vector amed vector of cell-> grp
#' @param newOrder vector of grp names in new order
#'
#' @return named vector new order
#' 
#' @export
reorderCellsByGrp<-function(
  nvect,# named vector of cell-> grp
  newOrder # vector of grp names in new order
){
  newAns<-vector()
  for(i in 1:length(newOrder)){
    nname<-newOrder[i]
    aib<-nvect[which(nvect==nname)]
    cat(nname," ", length(aib),"\n")
    newAns<-append(newAns, aib)
  }
  newAns
}

#' make tsne from pca
#'
#' make tsne from pca
#'
#' @param expDat expDat
#' @param recRes result of running gpaRecurse
#' @param perplexity (30)
#' @param theta (0.30)
#' @param weighted (TRUE) whether to use PCs from deeper lelves, and to weight them
#'
#' @return tsne matrix
#' 
#' @export
pca_to_tsne<-function
(expDat,
 recRes, # result of running gpaRecurse
 perplexity=30,
 theta=0.30,
 weighted=TRUE){

  tmpMat<-pca_project_all(expDat, recRes)
  if(weighted){
      
    weights<-rep(1, ncol(tmpMat))
    xi<- grep("L1_",colnames(tmpMat))
    weights[xi]<-10

    xi2<- grep("L2_",colnames(tmpMat))
    weights[xi2]<-.5

    xi3<- grep("L3_",colnames(tmpMat))
    weights[xi3]<-.25

    xi4<- grep("L4_",colnames(tmpMat))
    weights[xi4]<-.1
   
    datMat <- tmpMat %*% diag(weights)

  }
  else{
    ###datMat<-recRes$results[[1]]$gpRes$pcaRes$pcaRes$x
    datMat<-tmpMat
  }

  tres<-Rtsne(datMat, pca=FALSE, perplexity=perplexity, theta=theta, max_iter=2e3)
  xres<-tres$Y
  colnames(xres)<-c("TSNE.1", "TSNE.2")
  rownames(xres)<-rownames(datMat)
  xres
}


project_pca<-function
(expDat,
  pcRes,
  pcs=FALSE){

  if(length(pcs)==1){
    if(!pcs){
      pcs<-1:ncol(pcRes$rotation)
    }
  }
  scale(t(expDat), pcRes$center, pcRes$scale) %*% pcRes$rotation[,pcs]
}

pca_project_gpa<-function
(expDat,
 gpaRecRes,
 lName){

  myResult<-gpaRecRes$results[[lName]]$gpRes$pcaRes

  # get the varGenes
  varGenes<-myResult$varGenes

  # get the PCs
  pcs<-1:ncol(myResult$pcaRes$x)

  ans<-project_pca(expDat[varGenes,], myResult$pcaRes, pcs=pcs)
  ans
 }


pca_project_all<-function
(expDat,
gpaRecRes){
 
  resNames<-names(gpaRecRes$results)
  tmpAns<-list()
  colCount<-0

  allCnames<-vector()
  for(resName in resNames){
    cat(resName,"\n")
    tmpAns[[resName]]<-pca_project_gpa(expDat, gpaRecRes, resName)
    cnames<-paste(resName, "_",colnames(tmpAns[[resName]]), sep='')
    colnames(tmpAns[[resName]])<-cnames
    colCount<-colCount + ncol(tmpAns[[resName]])
    allCnames<-append(allCnames, cnames)
  }

  ans<-matrix(0, nrow=ncol(expDat), ncol=colCount)
  rownames(ans)<-colnames(expDat)
  stp<-0
  for(resName in resNames){
    str<-stp+1
    stp<-str + ncol(tmpAns[[resName]]) - 1
    ans[,str:stp]<-tmpAns[[resName]]
  }
  colnames(ans)<-allCnames
  ans

}











myGrpSort<-function(grps){
  grpNames<-unique(grps)
###   llevels<- as.numeric(unlist(lapply(strsplit(grpNames, "L"), "[[",2))) NEED TO FIX THIS TO repeat over Ls
  xix<-as.numeric(unlist(lapply(strsplit(grpNames, "_G"), "[[",2)))
   grpNames[order(xix)]
}

#' heatmap genes and groups
#'
#' heatmap genes and groups
#'
#' @param expDat expDat
#' @param genes genes
#' @param grps vector of cellnames -> grp label
#' @param maxPerGrp 100
#' @param cRow =FALSE,
#' @param cCol =FALSE,
#' @param limits =c(0,10),
#' @param toScale =FALSE,
#' @param fontsize_row =4
#'
#' @return pheatmap
#' 
#' @export
hm_gpa_sel<-function(
  expDat,
  genes,
  grps, ## vector of cellnames -> grp label
  maxPerGrp=100,
  cRow=FALSE,
  cCol=FALSE,
  limits=c(0,10),
  toScale=FALSE,
  fontsize_row=4,
  reOrderCells=FALSE){

  
  allgenes<-rownames(expDat)
  missingGenes<-setdiff(genes, allgenes)
  if(length(missingGenes)>0){
    cat("Missing genes: ", paste0(missingGenes, collapse=","), "\n")
    genes<-intersect(genes, allgenes)
  }

  value<-expDat[genes,]
  if(toScale){
      if(class(value)[1]!='matrix'){
        value <- t(scale(Matrix::t(value)))
      }
      else{
        value <- t(scale(t(value)))
      }
    }

  value[value < limits[1]] <- limits[1]
  value[value > limits[2]] <- limits[2]

  groupNames<-unique(grps)
  if(reOrderCells){
    grps<-grps[order(grps)]
    groupNames<-sort(unique(grps))
  }

  cells<-names(grps)  

##
 ## groupNames<-myGrpSort(grps)
##

  cells2<-vector()
  for(groupName in groupNames){
    xi<-which(grps==groupName)
    if(length(xi)>maxPerGrp){
      tmpCells<-sample(cells[xi], maxPerGrp)
    }
    else{
      tmpCells<-cells[xi]
    }
    cells2<-append(cells2, tmpCells)
  }
  value<-value[,cells2]

  xcol <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(groupNames))
    names(xcol) <- groupNames
    anno_colors <- list(group = xcol)

    xx<-data.frame(group=as.factor(grps))
    rownames(xx)<-cells

   val_col <- colorRampPalette(rev(brewer.pal(n = 12,name = "Spectral")))(25)
   #val_col <- colorRampPalette(brewer.pal(n = 12,name = "Spectral"))(100)

  pheatmap(value, cluster_rows = cRow, cluster_cols = cCol, color=val_col,
        show_colnames = FALSE, annotation_names_row = FALSE,
##        annotation_col = annTab,
        annotation_col = xx,
        annotation_names_col = FALSE, annotation_colors = anno_colors, fontsize_row=fontsize_row)
}




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
sc_hmClass<-function(
  classMat,
  grps, ## vector of cellnames -> grp label 
  isBig=FALSE,
  maxPerGrp=100,
  cRow=FALSE,
  cCol=FALSE,
  fontsize_row=4,
  scale=FALSE
){
 
  cools<-colorRampPalette(c("black", "limegreen", "yellow"))( 100 )
  bcol<-'white';
  if(isBig){
    bcol<-NA;
  }


  grps<-grps[order(grps)]
  cells<-names(grps)
  groupNames<-sort(unique(grps))

  cells2<-vector()
  for(groupName in groupNames){
    xi<-which(grps==groupName)
    if(length(xi)>maxPerGrp){
      tmpCells<-sample(cells[xi], maxPerGrp)
    }
    else{
      tmpCells<-cells[xi]
    }
    cells2<-append(cells2, tmpCells)
  }
  classMat<-classMat[,cells2]

  xcol <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(groupNames))
  names(xcol) <- groupNames
  anno_colors <- list(group = xcol)

  xx<-data.frame(group=as.factor(grps))
  rownames(xx)<-cells

  if(scale){
    mymin<-min(classMat)
    mymax<-max(classMat)
  }
  else{
    mymin<-0
    mymax<-1
  }
  
  pheatmap(classMat, col=cools, breaks=seq(from=mymin, to=mymax, length.out=100), cluster_rows = cRow, cluster_cols = cCol,
        show_colnames = FALSE, annotation_names_row = FALSE,
##        annotation_col = annTab,
        annotation_col = xx,
        annotation_names_col = FALSE, annotation_colors = anno_colors, fontsize_row=fontsize_row)

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
#' @param gpaRes gpRes
#' @param llevel llevel
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
#' @param gpaRes gpRes
#' @param legend whether to display it
#'
#' @return ggplot
#' 
#' @export
#'
plotGPApca<-function(xres, legend=FALSE)
{

  xdat<-xres$gpRes$pcaRes$pcaRes$x
  grps<-xres$bundleRes$result
  xdat<-xdat[names(grps),]
  aDat<-data.frame(pc1=xdat[,1], pc2=xdat[,2], group=as.character(grps) )
 
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
#' @param xtree xtree
#' @param legend whether to display it
#'
#' @return ggplot
#' 
#' @export
#'
plot_pca_gpa<-function(xtree, legend=FALSE)
{

  aDat<-data.frame(pc1=xtree$results[[1]]$gpRes$pcaRes$pcaRes$x[,1], pc2=xtree$results[[1]]$gpRes$pcaRes$pcaRes$x[,2], group=as.factor(xtree$groups)) 
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
#' @param gpaRes gpRES
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
plot_tsne<-function(sampTab, tsRes, cname="study_id", themeWhite=TRUE)
{
  xres<-cbind(sampTab, tsRes)
  xi<-which(colnames(xres)==cname)
  colnames(xres)[xi]<-"group"
  ColorRamp <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(unique(xres$group)))

  res<-ggplot(xres, aes(x=TSNE.1, y=TSNE.2, colour=group, fill=group) ) + geom_point(pch=21, alpha=3/4, size=1, stroke=0) + scale_colour_manual(values=ColorRamp) + scale_fill_manual(values=ColorRamp) #+ facet_wrap( ~ k, nrow=3)
  if(themeWhite){
    res <- res + theme_bw()
  }
  else{
    res <- res + theme_black()
  }
  res
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


#' plot tsne and genes
#'
#' plot tsne and genes 
#'
#' @param tsneDat tsne matrix
#' @param expDat expression matrix
#' @param genesToPlot genes to plot
#' @param colorPal "BuPu"
#' @param revCol TRUE
#' @param toScale scale?
#' @param limits
#' @return ggplot
#' 
#' @export
#'
tsneMultsimp<-function(
  tsneDat, # cols TSNE.1, TNSE.2
  expDat,
  genesToPlot, # genes to plot
  colorPal="BuPu",
  revCol=TRUE,
  toScale=TRUE,
  limits=c(-3,3),
  themeWhite=TRUE
 ){

  require(tidyr)

  allgenes<-rownames(expDat)

  missing<-setdiff(genesToPlot, allgenes)
  cat(paste0("missing ",missing,collapse=","),"\n")
  genesToPlot<-intersect(genesToPlot, allgenes)

  value<-expDat[genesToPlot,]

  if(toScale){
    value <- t(scale(t(value)))
  }
  value[value < limits[1]] <- limits[1]
  value[value > limits[2]] <- limits[2]


  tsne<-as.data.frame(tsneDat[,c("TSNE.1", "TSNE.2")])
  tsne<-cbind(tsne, t(value))

  tsneLong<-gather_(tsne, "gene", "expression", genesToPlot)
  tsneLong$gene_f = factor(tsneLong$gene, levels=genesToPlot)

  if(revCol){
    ColorRamp <- rev(colorRampPalette(rev(brewer.pal(n = 7,name = colorPal)))(100))[10:100]
  }
  else{
    ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7,name = colorPal)))(100)
  }
  res<- ggplot(tsneLong, aes(x=TSNE.1, y=TSNE.2, colour=expression) ) + 
    geom_point(pch=19, alpha=2/4, size=.25) + 
    scale_colour_gradientn(colours=ColorRamp) +
    facet_wrap( ~ gene_f)

  if(themeWhite){
    res<-res + theme_bw()
  }
  else{
    res<-res + theme_black()
 }
 res

#   tsneLong
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









hm_recRes<-function(
  expDat,
  recRes,
  toLevel=1,
  topx=10,
  maxPerGrp=100,
  toScale=FALSE,
  cRow=FALSE,
  cCol=FALSE,
  fontsize_row=4,
  limits=c(0,10))
{
  
  tGenes<-vector()
  nnames<-names(recRes$results)
  nodeNames<-vector()

  # find matching levels
  for(i in 1:toLevel){
    aa<-paste0("L",i)
    x<-nnames[grep(aa, nnames)]
    nodeNames<-append(nodeNames, x)
  }

  nodeNames<-sort(nodeNames)
  geneGrps<-vector()
  for(nodeName in nodeNames){
    #cat(nodeName,"\n")
    xres<-recRes$results[[nodeName]]
    diffExp<-xres$diffExp
    cluNames<-names(diffExp)
    cluNames<-sort(cluNames)
    ct1<-lapply( diffExp[cluNames], getTopGenes, topx)
    ct1<-unique(unlist(ct1))
    tGenes<-append(tGenes, ct1)
  }

  tGenes<-unique(tGenes)

  value<-expDat[tGenes,]
  if(toScale){
      value <- t(scale(t(value)))
  }

  value[value < limits[1]] <- limits[1]
  value[value > limits[2]] <- limits[2]

###  grps<-recRes$groups
  grps<-recRes$grp_list[[toLevel]]
  grps<-sort(grps)
  cells<-names(grps)
  groupNames<-unique(grps)

#  
  

###  bundleRes$result<-sort(bundleRes$result)

  cells2<-vector()
  for(groupName in groupNames){
    #cat(groupName,"\n")
    xi<-which(grps==groupName)
    if(length(xi)>maxPerGrp){
      tmpCells<-sample(cells[xi], maxPerGrp)
    }
    else{
      tmpCells<-cells[xi]
    }
    cells2<-append(cells2, tmpCells)
  }
  value<-value[,cells2]

  xcol <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(groupNames))
    names(xcol) <- groupNames
    anno_colors <- list(group = xcol)

    xx<-data.frame(group=as.factor(grps))
    rownames(xx)<-cells

##    yy<-data.frame(group=as.factor(geneGrps))
##    rownames(yy)<-tGenes


  pheatmap(value, cluster_rows = cRow, cluster_cols = cCol,
        show_colnames = FALSE, annotation_names_row = FALSE,
##        annotation_col = annTab,
        annotation_col = xx,
##        annotation_row = yy,
        annotation_names_col = FALSE, annotation_colors = anno_colors, fontsize_row=fontsize_row)


}

hm_gpa<-function(
  expDat,
  gpaRes,
  topx=10,
  maxPerGrp=100,
  toScale=FALSE,
  limits=c(0,10),
  fontsize_row=5)
{
  
  hm_diff(expDat, gpaRes$diffExp, gpaRes$bundleRes, topx=topx, maxPerGrp=maxPerGrp, toScale=toScale, limits=limits, fontsize_row=fontsize_row)
}


hm_diff<-function(
  expDat,
  diffRes,
  bundleRes,
  maxPerGrp=100,
  topx=10, 
  cRow=FALSE,
    cCol=FALSE,
    limits=c(0,10),
    toScale=FALSE,
  fontsize_row=3){
  
  ct1<-lapply( diffRes, getTopGenes, topx)
  ct1<-unique(unlist(ct1))
  value<-expDat[ct1,]
  if(toScale){
      value <- t(scale(t(value)))
    }

  value[value < limits[1]] <- limits[1]
    value[value > limits[2]] <- limits[2]

    bundleRes$result<-sort(bundleRes$result)


    cells<-names(bundleRes$result)

#   value<-value[,cells]
 
  groupNames<-unique(bundleRes$result)

  cells2<-vector()
  for(groupName in groupNames){
    cat(groupName,"\n")
    xi<-which(bundleRes$result==groupName)
    cat(length(xi),"\n")
    if(length(xi)>maxPerGrp){
      cat(maxPerGrp,"\n")
      tmpCells<-sample(cells[xi], maxPerGrp)
    }
    else{
      tmpCells<-cells[xi]
    }
    cells2<-append(cells2, tmpCells)
  }
  value<-value[,cells2]

  xcol <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(groupNames))
    names(xcol) <- groupNames
    anno_colors <- list(group = xcol)

    xx<-data.frame(group=as.factor(bundleRes$result))
    rownames(xx)<-cells

  pheatmap(value, cluster_rows = cRow, cluster_cols = cCol,
        show_colnames = FALSE, annotation_names_row = FALSE,
##        annotation_col = annTab,
        annotation_col = xx,
        annotation_names_col = FALSE, annotation_colors = anno_colors, fontsize_row=fontsize_row)
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

pcaPlot_recRes<-function(recRes, nodeName="L1_G1",legend=TRUE)
{

  gpaRes<-recRes$results[[nodeName]]$gpRes
  bundleRes<-recRes$results[[nodeName]]$bundleRes  
  aDat<-data.frame(pc1=gpaRes$pcaRes$pcaRes$x[,1], pc2=gpaRes$pcaRes$pcaRes$x[,2], group=as.character(bundleRes$result) )
  ColorRamp <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(unique(aDat$group)))
  if(legend){
    ans<-ggplot(aDat, aes(x=pc1, y=pc2, colour=group, pch=group) ) + geom_point(alpha=3/4, size=.5) + theme_bw() + scale_colour_manual(values=ColorRamp)
  }
  else{
    ans<-ggplot(aDat, aes(x=pc1, y=pc2, colour=group) ) + geom_point(pch=19, alpha=3/4, size=.5) + theme_bw() + scale_colour_manual(values=ColorRamp) + theme(legend.position="none")
  }
  ans
}

plot_bundle<-function(gpaRes, bundleRes, legend=TRUE)
{

  aDat<-data.frame(pc1=gpaRes$pcaRes$pcaRes$x[,1], pc2=gpaRes$pcaRes$pcaRes$x[,2], group=as.character(bundleRes$result) )
  ColorRamp <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(unique(aDat$group)))
  if(legend){
    ans<-ggplot(aDat, aes(x=pc1, y=pc2, colour=group) ) + geom_point(pch=19, alpha=3/4, size=.5) + theme_bw() + scale_colour_manual(values=ColorRamp)
  }
  else{
    ans<-ggplot(aDat, aes(x=pc1, y=pc2, colour=group) ) + geom_point(pch=19, alpha=3/4, size=.5) + theme_bw() + scale_colour_manual(values=ColorRamp) + theme(legend.position="none")
  }
  ans
}

# theme_dviz_hgrid  from https://github.com/clauswilke/dataviz/blob/master/R/themes.R
#' @export
theme_dviz_hgrid <- function(font_size = 14, font_family = "") { 
  color = "grey90"
  line_size = 0.5

  # Starts with theme_cowplot and then modify some parts
  theme_cowplot(font_size = font_size, font_family = font_family) %+replace% 
  theme( 
    # make horizontal grid lines
    panel.grid.major = element_line(colour = color, 
    size = line_size), 
    panel.grid.major.x = element_blank(), 
    
    # adjust axis tickmarks
    axis.ticks = element_line(colour = color, size = line_size), 

    # adjust x axis
    axis.line.x = element_line(colour = color, size = line_size), 

    # no y axis line
    axis.line.y = element_blank() ) 
} 



### theme_dviz_hgrid  FROM https://github.com/clauswilke/dataviz/blob/master/R/themes.R
theme_dviz_hgrid <- function(font_size = 14, font_family = "") { 
  color = "grey90"
  line_size = 0.5

# Starts with theme_cowplot and then modify some parts
  theme_cowplot(font_size = font_size, font_family = font_family) %+replace% 
  theme( 

    # make horizontal grid lines
    panel.grid.major = element_line(colour = color, 
    size = line_size), 
    panel.grid.major.x = element_blank(), 

    # adjust axis tickmarks
    axis.ticks = element_line(colour = color, size = line_size), 

    # adjust x axis
    axis.line.x = element_line(colour = color, size = line_size), 
    # no y axis line
    axis.line.y = element_blank() 
  ) 
}



# theme_black from https://gist.github.com/jslefche/eff85ef06b4705e6efbc
#' @export
theme_black = function(base_size = 12, base_family = "") { 
  library(gridExtra) 
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
  theme( 
    # Specify axis options
    axis.line = element_blank(), 
    axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9), 
    axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9), 
    axis.ticks = element_line(color = "white", size = 0.2), 
    #axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)), 
    #axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)), 
    axis.ticks.length = unit(0.3, "lines"), 
    # Specify legend options
    legend.background = element_rect(color = NA, fill = "black"), 
    legend.key = element_rect(color = "white", fill = "black"), 
    legend.key.size = unit(1.2, "lines"), 
    legend.key.height = NULL, 
    legend.key.width = NULL, 
    legend.text = element_text(size = base_size*0.8, color = "white"), 
    legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"), 
    legend.position = "right", 
    legend.text.align = NULL, 
    legend.title.align = NULL, 
    legend.direction = "vertical", 
    legend.box = NULL, 
    # Specify panel options
    panel.background = element_rect(fill = "black", color = NA), 
    panel.border = element_rect(fill = NA, color = "white"), 
    panel.grid.major = element_line(color = "grey35"), 
    panel.grid.minor = element_line(color = "grey20"), 
    #panel.margin = unit(0.5, "lines"), 
    # Specify facetting options
    ###strip.background = element_rect(fill = "grey30", color = "grey10"), 
    strip.background = element_blank(),
    ###strip.text.x = element_text(size = base_size*0.8, color = "white"), 
    strip.text.x = element_text(size = base_size, color = "white"), 
    strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90), 
    # Specify plot options
    plot.background = element_rect(color = "black", fill = "black"), 
    ###plot.title = element_text(size = base_size*1.2, color = "white"), 
    plot.title = element_text(size = base_size, color = "white"), 
    ##plot.margin = unit(rep(1, 4), "lines")
  ) 
}


