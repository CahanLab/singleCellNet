# Patrick Cahan (C) 2017
# patrick.cahan@gmail.com

# G: gene stats
# P: PCA
# A: assign cells to groups


snn_iter<-function
(knnRes,
### mydist,### will remove later, uselfu for silh now
 eps=seq(1,15,by=2),
 minPts=seq(5,30, by=5),
 k=30)
{
    combs<-expand.grid(eps,minPts)
    ansList<-list()
    sils<-vector()
    ks<-vector()
    for(i in seq(nrow(combs))){
        ansList[[i]]<-snn_clust(knnRes,k, eps=combs[i,1], minPts=combs[i,2])
        ks<-append(ks, length(unique(ansList[[i]])))
       ## sil<-silhouette(ansList[[i]], mydist)
       ## sils<-append(sils, summary(sil)$avg.width)
       ## cat("Eps: ",combs[i,1], " ... minPts: ",combs[i,2], " --> ",summary(sil)$avg.width, "\n")
    }
    ## params<-cbind(combs, silh=sils)
##     colnames(params)<-c("eps", "minPts", "Silh")
  
  multClusters<-which(ks!=1)

  override<-FALSE
  if(length(multClusters)==0){
    override<-TRUE
    params<-combs
  }
  else{
    ansList<-ansList[multClusters]
    names(ansList)<-1:length(ansList)
    params<-combs[multClusters,]
    params<-cbind(params,kval=ks[multClusters])
    colnames(params)<-c("eps", "minPts", "kval")
  }
  
  list(params=params, groups=ansList, override=override)
}



clearBorder<-function
(groups,
 knnObj,
 iters=10){

   for(i in 1:iters){
     xi<-which(groups==0)
     newGroups<-groups
   
    for(bPoint in xi){
        x<-groups[  knnObj$id[bPoint,] ]
        notZero<-which(x!=0)
        if(length(notZero)>0){
            newGroups[bPoint] <- as.numeric(names(which.max(table( x[notZero] ) ) ))
        }
      }
    groups<-newGroups
  }
  groups
}


snn_clust<-function
(knnRes,
 xDist,
 k=30,
 eps=1,
 minPts=15){
    tmpRes<-sNNclust(knnRes, k=k, eps=eps, minPts=minPts, borderPoints=TRUE)
    names(tmpRes$cluster)<-rownames(knnRes$id)
    ans<-clearBorder(tmpRes$cluster, knnRes)
    ans
}



bestGmc<-function(vect, grange){
  diffs<-diff(vect)
  grange[which.max(diffs)+1]
}

bestPC<-function(
  pcSDs,
  threshold=0.1){

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
  if(mName=="cutree"){
      ans<-A_cutree(gpRes, kvals=kvals)
  }
  if(mName=='mclust'){
    ans<-A_mclust(gpRes, kvals=kvals)
  }
  if(mName=='sNN_clust'){
    ans<-A_snnClust(gpRes, kNN_k=30, eps_range=seq(1,15,by=2),minPts_range=seq(5,30,by=5))
  }
  ans
}



# return a standard A_result (method=, kvals=, results=)
A_kmeans<-function(
	gpRes,
	kvals=2:5,
  max.iter=50,
  pam_pam=TRUE
){

	ans<-list()
	cellNames<-rownames(gpRes$pcaRes$pcaRes$x)
	for(kval in kvals){
		
   ###### 12-26-17 
   if(pam_pam){
      tmpRes<-pam(gpRes$pcaRes$pcaRes$x, kval, cluster.only=TRUE,pamonce=1)
		  x<-tmpRes
		  names(x)<-cellNames
    }
    else{
      tmpRes<-kmeans(gpRes$pcaRes$pcaRes$x, kval, iter.max=max.iter)
      x<-tmpRes$cluster
      names(x)<-cellNames
    }
		ans[[kval]]<-x
	}
	list(method="kmeans", kvals=kvals, results=ans, override=FALSE)
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

	list(method="cutree", kvals=kvals, results=ans, override=FALSE)
}



if(FALSE){
# don't select, let A_bundle select based on silhouette
A_mclust<-function(
	gpRes,
	kvals=2:5){

	ans<-list()

	for(kval in kvals){
		tmpRes<-Mclust(gpRes$pcaRes$pcaRes$x, G=kval)
		ans[[kval]]<-tmpRes$classification
	}
	list(method="Mclust", kvals=kvals, results=ans, override=FALSE)
}	
}


# Mclust select k based on BIC plateau
###12-27-17 if(FALSE){
A_mclust<-function(
  gpRes,
  ngenes=0,
  kvals=2:5
){

  cat("Mclust...\n")
  if(ngenes==0){
    ngenes<-nrow(gpRes$pcaRes$pcaRes$x)
  }
  ans<-list()
  bics<-mclustBIC(gpRes$pcaRes$pcaRes$x[1:ngenes,], G=kvals)#, modelNames=c("VVV"))
  myG<-bestGmc(bics[kvals-1],kvals) 
  tmpRes<-Mclust(gpRes$pcaRes$pcaRes$x[1:ngenes,], G=myG)#, modelNames=c("VVV"))
  ans[[myG]]<-tmpRes$classification
  cat(myG,"\n")
  list(method="Mclust", kvals=myG, results=ans, override=FALSE)
} 
###}

# mclust select based on MAC BIC
if(FALSE){
A_mclust<-function(
  gpRes,
  kvals=2:5){

  ans<-list()
  tmpRes<-Mclust(gpRes$pcaRes$pcaRes$x, G=kvals)
  kval<-tmpRes$G
  ans[[kval]]<-tmpRes$classification
  list(method="Mclust", kvals=kval, results=ans, override=FALSE)
} 
}


A_snnClust<-function(
  gpRes,
  kNN_k=30,
  eps_range=seq(1,15,by=2),
  minPts_range=seq(5,30,by=5)){

    aDist<-gpRes$xdist
    knnRes<-kNN(aDist, k=kNN_k)

##    tmpAns<-snn_iter(knnRes, aDist, eps=eps_range, minPts=minPts_range, k=kNN_k)

    tmpAns<-snn_iter(knnRes, eps=eps_range, minPts=minPts_range, k=kNN_k)

    ans<-list()
    kvals<-1:nrow(tmpAns[['params']])
    ans<-tmpAns$groups

    list(method="sNN_clust", kvals=kvals, results=ans, override=tmpAns$override)

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
  methods=c("mclust", "kmeans", "cutree", "sNN_clust"),
  handicap=0.00, # how much to reduce other methods sil in comparison to mclust
  adjust=TRUE # whether to penalize clustering with small clusters
){


  aRes_list<-list()
  silh_list<-list()
  method.winner<-''
  result.winner<-''
  k.winner<-''
  sil.clusters<-vector()

  maxes<-rep(0, length(methods))
  names(maxes)<-methods
  cat("----------------------------------------\n")
  for(mname in methods){
    
    if(mname=='sNN_clust'){
      tabs<-"\t"
    }
    else{
       tabs<-"\t\t"
    }

    aRes_list[[mname]] <- A_method(gpRes, kvals, mname)
    if(!aRes_list[[mname]]$override){
      silh_list[[mname]]<-A_silhouette(aRes_list[[mname]], gpRes$xdist, adjust=adjust)
      maxes[mname]<-max(silh_list[[mname]]$overall)
      cat(mname,tabs, maxes[mname],"\n")
    }
    else{
      silh_list[[mname]]<-list(overall = -1)
      maxes[mname]<-max(silh_list[[mname]]$overall)
      cat(mname,tabs, "-1\n")
    }

  }
  

  aai<-which(names(maxes)!='mclust')
  maxes[aai]<-maxes[aai]*(1-handicap)
  xi<-which.max(maxes)
  override<-FALSE
  if(maxes[xi]==0){
    override<-TRUE
    silhOA<-''
    silClusters<-''
    cat("Winner is none.\n")
  }
  else{
    method.winner<-methods[xi]

    sil.winner<-silh_list[[xi]]
    
    k.winner<-which.max(sil.winner$overall)

    silhOA<-sil.winner$overall[k.winner]
    silClusters<-sil.winner$cluster.sils[[k.winner]]
    

  #method.winner<-methods[xi]
    result.winner<-aRes_list[[xi]]$results[[k.winner]]
    winning_kval<-length(unique(result.winner))
    cat("Winner is ",method.winner," with ",winning_kval, " clusters.\n", sep='')
    cat("----------------------------------------\n")
  }

  list(override=override,method=method.winner, result=result.winner, k=k.winner, silh=list(overall=silhOA, sil.clusters=silClusters))

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

#  grpOrder<-myGrpSort(cellLabels)

#  specificSets[grpOrder]

  specificSets
}


geneEnr<-function(
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

  cat("Gene finding\n")
  ####myPatternG<-sc_sampR_to_pattern(as.vector(bundleRes$result))
  myPatternG<-sc_sampR_to_pattern(as.character(bundleRes$result))

  
  
 # aClust<-parallel::makeCluster(ncore, type='FORK')
  specificSets<-lapply(myPatternG, sc_testPattern, expDat=expDat)
 # stopCluster(aClust)
  
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
 methods=c("mclust", "kmeans", "cutree", "sNN_clust"),
 pcAuto=TRUE,
 adjust=TRUE){

  gpRes<-GP(expDat, nPCs=nPCs, dThresh=dThresh, zThresh=zThresh, meanType=meanType,  pcaMethod=pcaMethod, max.iter=max.iter, pcAuto=pcAuto)
	bundleRes<-A_bundle(gpRes, kvals=kvals, methods=methods, adjust=adjust)
  diffExp<-list()
  if(!bundleRes$override){
  ### bundleRes<-bundle_wrap(gpRes, kvals=kvals, methods=methods)
  	diffExp<-A_geneEnr(expDat[gpRes$pcaRes$varGenes,], gpRes, bundleRes, minSet=FALSE)
  }
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
  methods=c("mclust", "cutree", "kmeans", "sNN_clust"),
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


      tmpAns<-gpa(expAll[,cells_in_grp],
        kvals=kvals,
        nPCs=nPCs,
        dThresh=dThresh,
        zThresh=zThresh,
        meanType=meanType,
        max.iter=max.iter,
        methods=methods,
        pcaMethod=pcaMethod,
        adjust=adjust,
        pcAuto=pcAuto)


      # UPDATE notDone to done if 
      # if overall_silhs[[uniGrp]] > tmpAns$bundleRes$silh$overall
      # this indicates that subsequent clustering is not meaningful
###      if(tmpAns$bundleRes$silh$overall < min( overall_silhs[[uniGrp]], 0.5)){

  clusterFail<-TRUE
  if(!tmpAns$bundleRes$override){
    if( silMin ){
      clusterFail <- as.logical( min(tmpAns$bundleRes$silh$sil.clusters) <  (1-SilhDrop)*overall_silhs[[uniGrp]])
    }
    else{
        clusterFail <- as.logical(tmpAns$bundleRes$silh$overall <  (1-SilhDrop)*overall_silhs[[uniGrp]])
    }
  }
  if(clusterFail){
## 12-11-17-->    
      ##minCcount<-min(table(tmpAns$bundleRes$result))
       ## if(tmpAns$bundleRes$silh$overall <  (1-SilhDrop)*overall_silhs[[uniGrp]] || (minCcount < minClusterSize)){
       ### cat( "failed because ...\nmin: ",min(tmpAns$bundleRes$silh$sil.clusters),"\n")
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
#          cat(oldName," :: ",newName," :: ", length(xi),  "\n")
          cat("Cluster ",newName," has ", length(xi),  " cells.\n", sep='')
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
          ##	cat("bye bye to ... ",smallCluster,"\n")
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

