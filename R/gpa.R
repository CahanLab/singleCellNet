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
  	list(geneStats=geneStats, pcaRes=pcaRes)
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
	aDist,
	kvals=2:5
){
	atree<-hclust(aDist, "ward")
	ans<-list()
	for(kval in kvals){
		tmpRes<-cutree(atree, k=kval)
		ans[[kval]]<-tmpRes
	}

	list(method="cutree", kvals=kvals, results=ans)
}

A_mclust<-function(
	gpRes,
	kvals=2:5){

	ans<-list()
	for(kval in kvals){
		tmpRes<-Mclust(gpRes$pcaRes$pcaRes$x, G=kval)
		ans[[kval]]<-tmpRes$classification
	}
	list(method="Mclust", kvals=kvals, results=ans)
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
	kvals=2:5
){

	methods<-c("mclust", "kmeans", "cutree")
	xdist<-dist(gpRes$pcaRes$pcaRes$x)

##	cat("Mclust\n")
##	mres <-A_mclust(gpRes, kvals=kvals)
	cat("kmeans\n")
	kres <-A_kmeans(gpRes, kvals=kvals)
	cat("cutree\n")
	ctres<-A_cutree(xdist, kvals=kvals)


##	mcsilhs<-A_silhouette(mres, xdist)
	ksilhs <-A_silhouette(kres, xdist)
	ctsilhs<-A_silhouette(ctres, xdist)

##	mcmax<-max(mcsilhs$overall)
	kmax <-max(ksilhs$overall)
	ctmax<-max(ctsilhs$overall)

##	xi<-which.max(c(mcmax, kmax, ctmax))
##	sil.winner<-list(mcsilhs, ksilhs, ctsilhs)[[xi]]

	xi<-which.max(c( kmax, ctmax))
	sil.winner<-list(ksilhs, ctsilhs)[[xi]]
	k.winner<-which.max(sil.winner$overall)

	method.winner<-methods[xi]
##	result.winner<-list(mres, kres, ctres)[[xi]]$results[[k.winner]]
	result.winner<-list(kres, ctres)[[xi]]$results[[k.winner]]
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
 max.iter=20){


	gpRes<-GP(expDat, nPCs=nPCs, dThresh=dThresh, zThresh=zThresh, meanType=meanType,  pcaMethod=pcaMethod, max.iter=max.iter)
	bundleRes<-A_bundle(gpRes, kvals=kvals)
	diffExp<-A_geneEnr(expDat[gpRes$pcaRes$varGenes,], gpRes, bundleRes, minSet=FALSE)


if(FALSE){
  ans<-0
  ans<-simple_steam_cuttree(pcaRes$pcaRes$x[,1:nPCs],nClusters=nGrps)
  groups<-ans$grps
 # diffExp<-par_findSpecGenes(expDat[pcaRes$varGenes,], db_steamed$steamed$sampTab))
  pcs<-pcaRes$pcaRes$x[,1:nPCs]
  rownames(pcs)<-colnames(expDat)
  list(gs=geneStats, pcaRes=pcaRes, groups=groups, pcs=pcs, avesilh=ans$avesilh)
}
 
 list(gpRes=gpRes, bundleRes=bundleRes, diffExp=diffExp)

}


#' @export
gpa_recurse2<-function
(expAll,
  nPCs=2,
  sdExpl=0.01,
  dThresh=0,
  zThresh=2,
  meanType="overall_mean",
  nGrps=2,
  method="pcromp",
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
    tmpAns<-gpa_break(expAll, grps, isDone=isDone,nPCs=nPCs, minClusterSize=minClusterSize,sdExpl=sdExpl, dThresh=dThresh, zThresh=zThresh, meanType=meanType, nGrps=nGrps, method=method)
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
##  gpaRes,
  grps, # vector of cluster labels
  isDone, # bool vector indicating whether corresponding groups should be subjected to clustering
  minClusterSize=20,
  nPCs=2,
  dThresh=0,
  zThresh=2,
  meanType="overall_mean",
  kvals=2:5,
  pcaMethod="rpca",
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
      	pcaMethod=pcaMethod,
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

hm_gpa<-function(
	expDat,
	gpaRes,
	topx=10,
	maxPerGrp=100,
	toScale=FALSE,
	limits=c(0,10))
{
	
	hm_diff(expDat, gpaRes$diffExp, gpaRes$bundleRes, topx=topx, maxPerGrp=maxPerGrp, toScale=toScale, limits=limits)
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

