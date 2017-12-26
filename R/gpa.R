# Patrick Cahan (C) 2017
# patrick.cahan@gmail.com

# G: gene stats
# P: PCA
# A: assign cells to groups



#' project from pcs
#'
#' project from pcs
#'
#' @param recRes result of running gpaRecurse
#' @param perplexity (30)
#' @param theta (0.30)
#'
#' @return tsne matrix
#' 
#' @export



bestGmc<-function(vect, grange){
  diffs<-diff(vect)
  grange[which.max(diffs)+1]
}

bestPC<-function(
  pcSDs,
  threshold=0.05){

  ans<-min(which( abs(diff(pcSDs))<threshold ))
  if(ans==1){
    ans<-2
  }
  ans
  
}

GP<-function(
	expDat,
	dThresh=0, ### detection threshold for sc_statTab
	zThresh=2, ### zscore threshold for selecting variable genes
	meanType="overall_mean", ### basis to bin genes for variable gene detection
	pcaMethod="rpca", # rpca or prcomp
	nPCs=2, # top nPCs pcs to consider
  pcAuto=TRUE,
	max.iter=20
){
	cat("Calculating gene statistics...\n")
  	geneStats<-sc_statTab(expDat, dThresh=dThresh)
  	cat("PCA...\n")
  	pcaRes<-vg_pca(expDat, geneStats, zThresh=zThresh, meanType=meanType, method=pcaMethod, max.iter=max.iter)

    if(pcAuto){ # only valid if pcaMethod="prcomp"
      nPCs<-bestPC(pcaRes$pcaRes$sd)
      cat("Auto PCs: ",nPCs,"\n")
    }
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
	kvals=2:5,
  max.iter=50
){

	ans<-list()
	cellNames<-rownames(gpRes$pcaRes$pcaRes$x)
	for(kval in kvals){
		###tmpRes<-kmeans(gpRes$pcaRes$pcaRes$x, kval, iter.max=max.iter)
    tmpRes<-pam(gpRes$pcaRes$pcaRes$x, kval, cluster.only=TRUE,pamonce=1)
		x<-tmpRes
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


# don't select, let A_bundle select based on silhouette
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


# Mclust select k based on BIC plateau
if(FALSE){
A_mclust<-function(
  gpRes,
  ngenes=0,
  kvals=2:5
){

  if(ngenes==0){
    ngenes<-nrow(gpRes$pcaRes$pcaRes$x)
  }
  ans<-list()
  bics<-mclustBIC(gpRes$pcaRes$pcaRes$x[1:ngenes,], G=kvals)#, modelNames=c("VVV"))
  myG<-bestGmc(bics[kvals-1],kvals) 
  tmpRes<-Mclust(gpRes$pcaRes$pcaRes$x[1:ngenes,], G=myG)#, modelNames=c("VVV"))
  ans[[myG]]<-tmpRes$classification
  cat(myG,"\n")
  list(method="Mclust", kvals=myG, results=ans)
} 
}

# mclust select based on MAC BIC
if(FALSE){
A_mclust<-function(
  gpRes,
  kvals=2:5){

  ans<-list()
  tmpRes<-Mclust(gpRes$pcaRes$pcaRes$x, G=kvals)
  kval<-tmpRes$G
  ans[[kval]]<-tmpRes$classification
  list(method="Mclust", kvals=kval, results=ans)
} 
}

#return a list of kval-> {overall.silh=,__ , cluster.sils=__}
A_silhouette<-function(
	aRes,
	aDist,
  adjust=TRUE # whether to adjust overall silhouette based on clustering size

){

	kvals<-aRes$kvals
	overalls<-rep(0, max(kvals))
	indicies<-1:length(kvals)
	clusterList<-list()
  mins<-rep(0, max(kvals))
	for(i in seq(length(indicies))){

		kval <-kvals[i]
		sil<-silhouette(aRes$results[[kval]], aDist)
		overalls[kval]<-summary(sil)$avg.width
		clusterList[[kval]]<-summary(sil)$clus.avg.widths

    if(adjust){
      tmpSizes<-summary(sil)$clus.sizes
      adjTmp<-1+1/tmpSizes

      clusterList[[kval]]<-summary(sil)$clus.avg.widths / adjTmp

    }

    mins[kval]<-min(table(aRes$results[[kval]]))
	}
  if(adjust){
    adj<- 1 + 1/mins
    str<-min(kvals)
    stp<-max(kvals)
    overalls[str:stp]<-overalls[str:stp] / adj[str:stp]
  }
	list(overall=overalls, cluster.sils=clusterList)
}


# can be used to iterate over 2:nmax PCs
bundle_wrap<-function(
  gpRes,
   kvals=2:5,
  methods=c("mclust", "kmeans", "cutree"),
  handicap=0.00 # how much to reduce other methods sil in comparison to mclust
){

  
  gpX<-gpRes$pcaRes$pcaRes$x
  pcMax<-ncol(gpX)
  bList<-list()
  bScore<-vector()
  myrange<-2:pcMax
  for(stp in  myrange){
    gpTmp<-gpRes
    gpTmp$pcaRes$pcaRes$x<-gpX[,1:stp]
    gpTmp$xdist<-dist(gpX[,1:stp])
    
    tmpRes<-A_bundle(gpTmp, kvals=kvals, methods=methods, handicap=handicap)
    silO<-tmpRes$silh$overall
    bScore<-append(bScore, silO)
    bList[[stp]]<-tmpRes
    cat("testing pcas: ",stp,": ", silO ,"------\n")
  }
  cat("Winnder pcas: ",max(bScore)," ", myrange[which.max(bScore)],"\n")
  bList[[ myrange[which.max(bScore)] ]]

}

A_bundle<-function(
  gpRes,
  kvals=2:5,
  methods=c("mclust", "kmeans", "cutree"),
  handicap=0.00, # how much to reduce other methods sil in comparison to mclust
  adjust=TRUE # whether to penalize clustering with small clusters
){


  aRes_list<-list()
  silh_list<-list()
  maxes<-rep(0, length(methods))
  names(maxes)<-methods
  for(mname in methods){
    
    aRes_list[[mname]] <- A_method(gpRes, kvals, mname)
    silh_list[[mname]]<-A_silhouette(aRes_list[[mname]], gpRes$xdist, adjust=adjust)


    maxes[mname]<-max(silh_list[[mname]]$overall)
    cat(mname," (",which.max(silh_list[[mname]]$overall),")\t", maxes[mname],"\n")
  }

  aai<-which(names(maxes)!='mclust')
  maxes[aai]<-maxes[aai]*(1-handicap)
  xi<-which.max(maxes)
  method.winner<-methods[xi]

  sil.winner<-silh_list[[xi]]
  k.winner<-which.max(sil.winner$overall)
  cat("Winner:: ",k.winner," ",method.winner,"\n")

  #method.winner<-methods[xi]
  result.winner<-aRes_list[[xi]]$results[[k.winner]]

  list(method=method.winner, result=result.winner, k=k.winner, silh=list(overall=sil.winner$overall[k.winner], sil.clusters=sil.winner$cluster.sils[[k.winner]]))

}


#' find genes higher in a cluster compared to all other cells
#'
#' ind genes higher in a cluster compared to all other cells
#'
#' @param expDat expDat
#' @param cellLabels named vector of cell groups
#'
#' @return list of diffExp data framnes
#' 
#' @export
gnrAll<-function(
  expDat,
  cellLabels){

  myPatternG<-sc_sampR_to_pattern(as.character(cellLabels))
  specificSets<-lapply(myPatternG, sc_testPattern, expDat=expDat)
  cat("Done testing\n")

  grpOrder<-myGrpSort(cellLabels)

  specificSets[grpOrder]
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
  ####myPatternG<-sc_sampR_to_pattern(as.vector(bundleRes$result))
  myPatternG<-sc_sampR_to_pattern(as.character(bundleRes$result))

  
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
 pcaMethod="prcomp",
 max.iter=30,
 methods=c("mclust", "kmeans", "cutree"),
 pcAuto=TRUE,
 adjust=TRUE){

cat("nPCs ",nPCs,"\n") 
  gpRes<-GP(expDat, nPCs=nPCs, dThresh=dThresh, zThresh=zThresh, meanType=meanType,  pcaMethod=pcaMethod, max.iter=max.iter, pcAuto=pcAuto)
	bundleRes<-A_bundle(gpRes, kvals=kvals, methods=methods, adjust=adjust)
  ### bundleRes<-bundle_wrap(gpRes, kvals=kvals, methods=methods)
	diffExp<-A_geneEnr(expDat[gpRes$pcaRes$varGenes,], gpRes, bundleRes, minSet=FALSE)
  #names(diffExp)<-unique(bundleRes$result)
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


# run gpa while notDone==TRUE
# give clusters good names
# store results in a tree

#' find clusters, agnostically (kind of)
#'
#'find clusters, agnostically (kind of). Run PAM and cutree, across range of ks, select best based on silhouette
#'
#' @param expAll expDat
#' @param kvals kvals 2:5
#' @param nPCs number of PCs to use (only used if pcAuto=FALSE)
#' @param dThresh 0,
#' @param zThresh 2,
#' @param meanType ="overall_mean",
#' @param maxLevel maxLevel =4,
#' @param max.iter 30,
#' @param SilhDrop 0.25,
#' @param minClusterSize =42,
#' @param methods  methods =c("mclust", "cutree", "kmeans"),
#' @param pcaMethod  "prcomp",
#' @param numGenes (5)
#' @param silMin terminate clustering on a set of cells if the current min silh is < threshold. in contrast to using mean silh
#' @param pcAuto TRUE auto select the best PCs based on difference in Stdev 
#' @param adjust adjust TRUE whether to adjust silhouette based on the cluster size
#'
#' @return a very complex object
#' 
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
  pcaMethod="prcomp",
  numGenes=5,
  silMin=TRUE,
  pcAuto=TRUE,
  adjust=TRUE){
  

  set.seed(42)

  topNodeName<-"L1_G1"
  silhs<-list()
  silhs[[topNodeName]]<- -1

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

  overall_silhs[[topNodeName]]<- -1

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


      cat("GROUP!!!! ->>>> ",uniGrp,"\n")
      tmpAns<-gpa(expAll[,cells_in_grp],
        kvals=kvals,
        nPCs=nPCs,
        dThresh=dThresh,
        zThresh=zThresh,
        meanType=meanType,
        max.iter=max.iter,
        methods=methods,
        pcaMethod=pcaMethod,
        adjust=adjust)


      # UPDATE notDone to done if 
      # if overall_silhs[[uniGrp]] > tmpAns$bundleRes$silh$overall
      # this indicates that subsequent clustering is not meaningful
###      if(tmpAns$bundleRes$silh$overall < min( overall_silhs[[uniGrp]], 0.5)){

  clusterFail<-TRUE
  if( silMin ){
    clusterFail <- as.logical( min(tmpAns$bundleRes$silh$sil.clusters) <  (1-SilhDrop)*overall_silhs[[uniGrp]])
  }
  else{
      clusterFail <- as.logical(tmpAns$bundleRes$silh$overall <  (1-SilhDrop)*overall_silhs[[uniGrp]])
  }
    if(clusterFail){
## 12-11-17-->    
      ##minCcount<-min(table(tmpAns$bundleRes$result))
       ## if(tmpAns$bundleRes$silh$overall <  (1-SilhDrop)*overall_silhs[[uniGrp]] || (minCcount < minClusterSize)){
        cat( "failed because ...\nmin: ",min(tmpAns$bundleRes$silh$sil.clusters),"\n")
        notDone[cells_in_grp] <- 0
        grps[cells_in_grp] <- uniGrp
        ansList[[uniGrp]]<-NULL
      }
      else{
      # update group names
      # update silhouette
        tmpGrps<-tmpAns$bundleRes$result
        ####  ##oldNames<-unique(tmpGrps)
        oldNames<-names(tmpAns$diffExp)
        #tmpAns$diffExp)<-oldNames
        nnames<-vector()
        for(oldName in oldNames){
          newName<-paste0("L",count_level,"_G",group_count)
          nnames<-append(nnames, newName)
          group_count<-group_count+1
          xi<-names(which(tmpGrps==oldName))
          cat(oldName," :: ",newName," :: ", length(xi),  "\n")
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
        cat("testing here for group length hyp ...\n")
        grpLengths<-table( tmpGrps )
        if( any(grpLengths < minClusterSize )){
          smallClusters<-names(which(grpLengths<minClusterSize))
          for(smallCluster in smallClusters){
          	cat("bye bye to ... ",smallCluster,"\n")
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

