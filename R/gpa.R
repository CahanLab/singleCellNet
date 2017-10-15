# Patrick Cahan (C) 2017
# patrick.cahan@gmail.com

# G: gene stats
# P: PCA
# A: assign cells to groups

GP<-function(
	expDat,
	dThresh=0, ### detection threshold for sc_statTab
	zThresh=2, ### zscore threshold for selecting variable genes
	meanType="overall_mean", ### basis to bin genes for variable gene detection
	pcaMethod="rpca", # rpca or prcomp
	nPCs=2, # top nPCs pcs to consider
	max.iter=20
){
	cat("Calculating gene statistics...\n")
  	geneStats<-sc_statTab(expDat, dThresh=dThresh)
  	cat("PCA...\n")
  	pcaRes<-vg_pca(expDat, geneStats, zThresh=zThresh, meanType=meanType, method=pcaMethod, max.iter=max.iter)
  	pcaRes$pcaRes$x<-pcaRes$pcaRes$x[,1:nPCs]

    xdist<-dist(pcaRes$pcaRes$x)
  	list(geneStats=geneStats, pcaRes=pcaRes, xdist=xdist)
 }


# 
A_method<-function(
   gpRes,
   kvals=2:5,
   mName)
{
  if(mName=="kmeans"){
    ans<-A_kmeans(gpRes, kvals=kvals)
  }
  else{
    if(mName=="cutree"){
      ans<-A_cutree(gpRes, kvals=kvals)
    }
    else{
        ans<-A_mclust(gpRes, kvals=kvals)
      }
  }
  ans
}



# return a standard A_result (method=, kvals=, results=)
A_kmeans<-function(
	gpRes,
	kvals=2:5
){

	ans<-list()
	cellNames<-rownames(gpRes$pcaRes$pcaRes$x)
	for(kval in kvals){
		tmpRes<-kmeans(gpRes$pcaRes$pcaRes$x, kval)
		x<-tmpRes$cluster
		names(x)<-cellNames
		ans[[kval]]<-x
	}
	list(method="kmeans", kvals=kvals, results=ans)
}


A_cutree<-function(
	gpRes,
	kvals=2:5,
  linkage="ward.D"
#  linkage="average"
){
	atree<-hclust(gpRes$xdist, linkage)
	ans<-list()
	for(kval in kvals){
		tmpRes<-cutree(atree, k=kval)
		ans[[kval]]<-tmpRes
	}

	list(method="cutree", kvals=kvals, results=ans)
}

if(FALSE){A_mclust<-function(
	gpRes,
	kvals=2:5){

	ans<-list()
	for(kval in kvals){
		tmpRes<-Mclust(gpRes$pcaRes$pcaRes$x, G=kval)
		ans[[kval]]<-tmpRes$classification
	}
	list(method="Mclust", kvals=kvals, results=ans)
}	
}

A_mclust<-function(
  gpRes,
  kvals=2:5){

  ans<-list()
  tmpRes<-Mclust(gpRes$pcaRes$pcaRes$x, G=kvals)
  kval<-tmpRes$G
  ans[[kval]]<-tmpRes$classification
  list(method="Mclust", kvals=kval, results=ans)
} 


#return a list of kval-> {overall.silh=,__ , cluster.sils=__}
A_silhouette<-function(
	aRes,
	aDist
){

	kvals<-aRes$kvals
	overalls<-rep(0, max(kvals))
	indicies<-1:length(kvals)
	clusterList<-list()
	for(i in seq(length(indicies))){

		kval <-kvals[i]
		sil<-silhouette(aRes$results[[kval]], aDist)
		overalls[kval]<-summary(sil)$avg.width
		clusterList[[kval]]<-summary(sil)$clus.avg.widths
	}
	list(overall=overalls, cluster.sils=clusterList)
}

A_bundle<-function(
	gpRes,
	kvals=2:5,
  methods=c("mclust", "kmeans", "cutree")
){


  aRes_list<-list()
  silh_list<-list()
  maxes<-rep(0, length(methods))
  names(maxes)<-methods
  for(mname in methods){
    cat(mname,"\n")
    aRes_list[[mname]] <- A_method(gpRes, kvals, mname)
    silh_list[[mname]]<-A_silhouette(aRes_list[[mname]], gpRes$xdist)
    maxes[mname]<-max(silh_list[[mname]]$overall)
  }



  if(FALSE){
  	cat("Mclust\n")
  	mres <-A_mclust(gpRes, kvals=kvals)
  	cat("kmeans\n")
  	kres <-A_kmeans(gpRes, kvals=kvals)
  	cat("cutree\n")
  	ctres<-A_cutree(xdist, kvals=kvals)

  	mcmax<-max(mcsilhs$overall)
  	kmax <-max(ksilhs$overall)
  	ctmax<-max(ctsilhs$overall)

    cat("MClust max: ", mcmax,"\n")
    cat("kmeans max: ", kmax,"\n")
    cat("cutree max: ", ctmax,"\n")
  }

	xi<-which.max(maxes)
	method.winner<-methods[xi]


	sil.winner<-silh_list[[xi]]
	k.winner<-which.max(sil.winner$overall)
  cat("Winner:: ",k.winner," ",method.winner,"\n")

	method.winner<-methods[xi]
	result.winner<-aRes_list[[xi]]$results[[k.winner]]

	list(method=method.winner, result=result.winner, k=k.winner, silh=list(overall=sil.winner$overall[k.winner], sil.clusters=sil.winner$cluster.sils[[k.winner]]))
}

A_geneEnr<-function(
	expDat,
	gpaRes, 
	bundleRes,
	prune=FALSE, ### limit to genes exclusively detected as CT in one CT
  	minSet=TRUE,
  	ncore=4
){


if(FALSE){
  newST<-sampTab
  if(minSet){
    cat("Making reduced sampleTable\n")
    newST<-minTab(sampTab, dLevel)
  }

}  
  ans<-list()

  expDat<-expDat[, names(bundleRes$result)]

  cat("Making patterns\n")
  myPatternG<-sc_sampR_to_pattern(as.vector(bundleRes$result))

  
  cat("Testing patterns\n")
 # aClust<-parallel::makeCluster(ncore, type='FORK')
  specificSets<-lapply(myPatternG, sc_testPattern, expDat=expDat)
 # stopCluster(aClust)
  cat("Done testing\n")
  specificSets
}


#' @export
gpa<-function(expDat,
 kvals=2:5,
 nPCs=2,
 dThresh=0,
 zThresh=2,
 meanType="overall_mean",
 pcaMethod="rpca",
 max.iter=30,
 methods=c("mclust", "kmeans", "cutree")){

  gpRes<-GP(expDat, nPCs=nPCs, dThresh=dThresh, zThresh=zThresh, meanType=meanType,  pcaMethod=pcaMethod, max.iter=max.iter)
	bundleRes<-A_bundle(gpRes, kvals=kvals, methods=methods)
	diffExp<-A_geneEnr(expDat[gpRes$pcaRes$varGenes,], gpRes, bundleRes, minSet=FALSE)
  list(gpRes=gpRes, bundleRes=bundleRes, diffExp=diffExp)

}


makeNode<-function(
  nodeName,
  ncells,
  silh,
  topGenes)
{
  tmpNode<-Node$new(nodeName)
  tmpNode$cells<-ncells
  tmpNode$silh<-silh
  tmpNode$topGenes<-paste(topGenes, collapse=", ")
  tmpNode
}


# run gpa_break while notDone==TRUE
# give clusters good names
# store results in a tree
#' @export
gpa_recurse<-function(
  expAll,
  kvals=2:5,
  nPCs=2,
  dThresh=0,
  zThresh=2,
  meanType="overall_mean",
  maxLevel=4,
  max.iter=30,
  SilhDrop=0.25,
  minClusterSize=42,
  methods=c("mclust", "cutree", "kmeans"),
  numGenes=5){
  

  topNodeName<-"L1_G1"
  silhs<-list()
  silhs[[topNodeName]]<-0

  myTree<-makeNode(topNodeName, ncol(expAll), silh=0, topGenes="")
  
  listOfNodes<-list()
  listOfNodes[[topNodeName]]<-myTree


  grps<-rep(topNodeName, ncol(expAll))
  notDone<-rep(1, ncol(expAll))

  names(grps)<-colnames(expAll)
  names(notDone)<-colnames(expAll)

  ansList<-list()
  grp_list<-list()
  overall_silhs<-list()

  overall_silhs[[topNodeName]]<-0

  grp_list[[1]]<-grps

  gene_list<-list()

  count_level <- 2

  orderedGrps<-vector()

  while( (count_level <= maxLevel) && any(as.logical(notDone))){
    cat("Level: ",count_level,"\n")

    cellsToCluster<-names(which(notDone==TRUE))
    grpsToCluster<-grps[cellsToCluster]

    uniGrps<-unique(grpsToCluster)

    group_count<-1

    for(i in 1:length(uniGrps)){
      uniGrp<-uniGrps[i]

      currNode<-listOfNodes[[uniGrp]]

      cells_in_grp<-names(which(grps==uniGrp))
      tmpAns<-gpa(expAll[,cells_in_grp],
        kvals=kvals,
        nPCs=nPCs,
        dThresh=dThresh,
        zThresh=zThresh,
        meanType=meanType,
        max.iter=max.iter,
        methods=methods)


      # UPDATE notDone to done if 
      # if overall_silhs[[uniGrp]] > tmpAns$bundleRes$silh$overall
      # this indicates that subsequent clustering is not meaningful
     ### if(tmpAns$bundleRes$silh$overall < min( overall_silhs[[uniGrp]], 0.5)){
      if(tmpAns$bundleRes$silh$overall <  (1-SilhDrop)*overall_silhs[[uniGrp]]){
        notDone[cells_in_grp] <- 0
        grps[cells_in_grp] <- uniGrp
        ansList[[uniGrp]]<-NULL
      }
      else{
      # update group names
      # update silhouette
        tmpGrps<-tmpAns$bundleRes$result
        oldNames<-unique(tmpGrps)
        nnames<-vector()
        for(oldName in oldNames){
          newName<-paste0("L",count_level,"_G",group_count)
          nnames<-append(nnames, newName)
          cat(newName,"\n")
          group_count<-group_count+1
          xi<-names(which(tmpGrps==oldName))
          tmpGrps[xi]<-newName
          overall_silhs[[newName]] <- tmpAns$bundleRes$silh$sil.clusters[ oldName ]

          # create new node and add to tree ...
          # ... add data.tree code here ...
          xy_genes<-getTopGenes(tmpAns$diffExp[[oldName]],numGenes)
          ###newNode<-makeNode(newName, length(xi), silh=overall_silhs[[newName]], topGenes=getTopGenes(tmpAns$diffExp[[oldName]], 3))
          newNode<-makeNode(newName, length(xi), silh=overall_silhs[[newName]], topGenes=xy_genes)
          newNode<-currNode$AddChildNode(newNode)
          listOfNodes[[newName]]<-newNode
          gene_list[[newName]]<-xy_genes
        } 
        names(tmpAns$diffExp)<-nnames
        tmpAns$bundleRes$result<-tmpGrps

        # UPDATE notDone those cells in clusters with sizes < minClusterSize
        grpLengths<-table( tmpGrps )
        if( any(grpLengths < minClusterSize )){
          smallClusters<-names(which(grpLengths<minClusterSize))
          for(smallCluster in smallClusters){
            xi<-names(which(tmpGrps==smallCluster))
            notDone[xi]<-0
          }
        }
        grps[cells_in_grp]<-tmpGrps[cells_in_grp]
        ansList[[uniGrp]]<-tmpAns
      }      
    }
    grp_list[[count_level]]<-grps
    count_level<-count_level + 1
  }
  list(results=ansList, groups=grps, grp_list=grp_list, groupTree=myTree, geneList=gene_list)
}


#' @export
gpa_recurse2<-function(
  expAll,
  nPCs=2,
  sdExpl=0.01,
  dThresh=0,
  zThresh=2,
  meanType="overall_mean",
  nGrps=2,
  max=10,
  minClusterSize=42,
  ariThresh=.95){

  
  grps<-rep("0", ncol(expAll))
  isDone<-rep(0, ncol(expAll))

  ansList<-list()
  grp_list<-list()

  count_i<-1
  while(count_i <= max){
    cat(count_i,"\n")
    tmpAns<-gpa_break(expAll, 
      grps, 
      isDone=isDone,
      nPCs=nPCs, 
      minClusterSize=minClusterSize,
      dThresh=dThresh,
      zThresh=zThresh,
      meanType=meanType,
      nGrps=nGrps)
    xgrps<-tmpAns$groups
    grp_list[[count_i]]<-xgrps
    ari<-adjustedRandIndex(grps, xgrps)
    cat("Round ", count_i," ARI = ", ari,"\n")
    if(ari>ariThresh | all(isDone==TRUE) ){
      break
    }
    else{
      ansList[[count_i]]<-tmpAns
      grps<-xgrps
      isDone<-tmpAns$isDone
      count_i<-count_i+1
    }
  }
  list(results=ansList, groups=grps, grp_list=grp_list)
}



#' @export
# run gpa on each grp
gpa_break<-function(
  expDat,
  grps, # vector of cluster labels
  isDone, # bool vector indicating whether corresponding groups should be subjected to clustering
  minClusterSize=20,
  nPCs=2,
  dThresh=0,
  zThresh=2,
  meanType="overall_mean",
  kvals=2:5,
  max.iter=20
 ){


  ans_grp<-grps
  ans_list<-list()
  uniGrps<-unique(grps)
  ##dontCluster<-unique(grps[which(isDone==TRUE)])
  ##uniqGrps<-setdiff(uniGrps, dontCluster)

##  pcs<-matrix(0, nrow=ncol(expDat), ncol=2)
 ## rownames(pcs)<-colnames(expDat)

  for(i in 1:length(uniGrps)){
    uniGrp<-uniGrps[i]
    xi<-which(grps==uniGrp)
   ## if(length(xi)>=minClusterSize){      
      cat(uniGrp,":" ,length(xi),"\n")
    ###ans_list[[uniGrp]]<-gpa(expDat[,xi], minClusterSize=minClusterSize,nPCs=nPCs,sdExpl=sdExpl,dThresh=dThresh,zThresh=zThresh,meanType=meanType,gMax=gMax)
     
     ##**## ans_list[[uniGrp]]<-gpa_tree(expDat[,xi],nPCs=nPCs,sdExpl=sdExpl,dThresh=dThresh,zThresh=zThresh,meanType=meanType,nGrps=nGrps, method=method)

      gpaRes<-gpa(expDat[,xi],
      	nPCs=nPCs,
      	dThresh=dThresh,
      	zThresh=zThresh,
      	meanType=meanType,kvals=kvals,
      	max.iter=max.iter)
      
      	ans_list[[i]]<-gpaRes
if(FALSE){
      if( (any(table(gpaRes$groups) < minClusterSize))){ ###|| ( gpaRes$avesilh < aveSilhs[[uniqGrp]])){ 
        ans_list[[uniGrp]]<-""
        ans_grp[xi]<-uniGrp
        isDone[xi]<-1
      }
      else{
        snames<-colnames(expDat[,xi])
        pcs[snames,1:nPCs]<-gpaRes$pcs[,1:nPCs]
        ans_list[[uniGrp]]<-gpaRes
        ans_grp[xi]<-paste0(uniGrp,"_",ans_list[[uniGrp]]$groups)
      }
    }
}
#  list(groups=ans_grp, res=ans_list, isDone=isDone, pcs=pcs)#, aveSilhs=aveSils)
ans_list
}



######################
# plotting functions #
######################


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

#  	value<-value[,cells]
 
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

