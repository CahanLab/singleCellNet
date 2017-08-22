# Patrick Cahan (C) 2017
# patrick.cahan@gmail.com

# functions to assign samples/cells to groups using dbscan, mclust, treecut and frienemies

#' @export
cutTreeToNames<-function # converts the result of running cutreeDynamicTree to a list of entities
(entNames, # vector of names in same order as ctRes
 ctRes # result of running cutreeDynamicTree 
){
  ans<-list()
  grpNames<-unique(ctRes)
  for(grpName in grpNames){
    ans[[as.character(grpName)]]<-sort(entNames[ which(ctRes==grpName) ])
  }
  ans
}

if(FALSE){
#' @export
steam_treecut<-function
(sampTab,
 datMat,
 distanceMetric="spearman",
 linkageMethod="complete",
 nClusters=c(5,20),
 sClusters=c(10,NA),
 mmsizeRange=seq(25,500,by=25),
 deepSplitRange=0:4
 ){

	cat("Computing cell-to-cell distances\n")
	if(distanceMetric=='euclidean'){
		cors<-dist(datMat) 
	}
	else{
		cors<-cor(t(datMat), method=distanceMetric)
		cors<-as.dist(1-cors)
	}

	cat("Clustering\n")
	hcRes<-hclust(cors, linkageMethod)

 	if(is.na(sClusters[2])){
 		sClusters[2]<-ceiling(nrow(sampTab)*.75)
 	}
 	cat(sClusters[1],"<->",sClusters[2],"\n")

 	args<-expand.grid(mmsizeRange, deepSplitRange)
    nClus<-rep(0, nrow(args))
    minSizes<-rep(0, nrow(args))
    maxSizes<-rep(0, nrow(args))
 	for(i in seq(nrow(args))){
 		tmpAns<-.steam_dt(sampTab, hcRes, cors,minModSize=args[i,1],deepSplit=args[i,2])
 		grps<-table(tmpAns$group)
 		nClus[i]<-length(grps)
 		minSizes[i]<-min(grps)
 		maxSizes[i]<-max(grps)
 	}
 	ans<-data.frame(num_groups=nClus, smallest=minSizes, largest=maxSizes)
 	tmpAns<-cbind(args, ans)
# 	rownames(tmpAns)<-as.character(1:nrow(tmpAns))

 	# select the ones that fit criteria
 	#crit1<-rownames(tmpAns[tmpAns$num_groups> nClusters[1] & tmpAns$num_groups< nClusters[2],])
 	crit1<-which(tmpAns$num_groups> nClusters[1] & tmpAns$num_groups< nClusters[2])
 	crit2<-which(tmpAns$smallest>= sClusters[1] & tmpAns$largest<= sClusters[2])
 	x1<-intersect(crit1, crit2)
 
 	if(length(x1)>=1){
 		ai<-which.min(tmpAns[x1,]$num_groups)
 		#x<-which.min(tmpAns[x1[ai],]$largets)

 		x<-x1[ai]
 		if(length(x)>1){
 			x<-x[1]
 		}
 		ans<-.steam_dt(sampTab, hcRes, cors, minModSize=args[x,1], deepSplit=args[x,2])
 		cat("minMidSize=",args[x,1], "deepSplit=",args[x,2],"\n")
 		opt_params<-list(minMidSize=args[x,1], deepSplit=args[x,2])
 	}
 	else{
 		ans<-NULL
 		cat("Nothing fits criteria\n")
 		opt_params<-list(minMidSize=NA, deepSplit=NA)
 	}
 	ans<-ans[order(ans$group),]
 	args<-as.list(match.call())
 	
 	list(sampTab=ans, args=args, opt_params=opt_params)
}
}

#' @export
steam_treecut<-function
(sampTab,
 datMat,
 distanceMetric="spearman",
 linkageMethod="complete",
 nClusters=c(5,20),
 sClusters=c(10,NA),
 mmsizeRange=seq(25,500,by=25),
 deepSplitRange=0:4
 ){

	cat("Computing cell-to-cell distances\n")
	if(distanceMetric=='euclidean'){
		cors<-dist(datMat) 
	}
	else{
		cors<-cor(t(datMat), method=distanceMetric)
		cors<-as.dist(1-cors)
	}

	cat("Clustering\n")
	hcRes<-hclust(cors, linkageMethod)

 	if(is.na(sClusters[2])){
 		sClusters[2]<-ceiling(nrow(sampTab)*.75)
 	}
 	cat(sClusters[1],"<->",sClusters[2],"\n")

 	args<-expand.grid(mmsizeRange, deepSplitRange)
    nClus<-rep(0, nrow(args))
    minSizes<-rep(0, nrow(args))
    maxSizes<-rep(0, nrow(args))
 	for(i in seq(nrow(args))){
 		tmpAns<-.steam_dt(sampTab, hcRes, cors,minModSize=args[i,1],deepSplit=args[i,2])
 		grps<-table(tmpAns$group)
 		nClus[i]<-length(grps)
 		minSizes[i]<-min(grps)
 		maxSizes[i]<-max(grps)
 	}
 	ans<-data.frame(num_groups=nClus, smallest=minSizes, largest=maxSizes)
 	tmpAns<-cbind(args, ans)
# 	rownames(tmpAns)<-as.character(1:nrow(tmpAns))


 	# select the ones that fit criteria
 	#crit1<-rownames(tmpAns[tmpAns$num_groups> nClusters[1] & tmpAns$num_groups< nClusters[2],])
 	crit1<-which(tmpAns$num_groups> nClusters[1] & tmpAns$num_groups< nClusters[2])
 	crit2<-which(tmpAns$smallest>= sClusters[1] & tmpAns$largest<= sClusters[2])
 	x1<-intersect(crit1, crit2)
 
 	if(length(x1)>=1){
 		ai<-which.min(tmpAns[x1,]$num_groups)
 		#x<-which.min(tmpAns[x1[ai],]$largets)

 		x<-x1[ai]
 		if(length(x)>1){
 			x<-x[1]
 		}
 		ans<-.steam_dt(sampTab, hcRes, cors, minModSize=args[x,1], deepSplit=args[x,2])
 		cat("minMidSize=",args[x,1], "deepSplit=",args[x,2],"\n")
 		opt_params<-list(minMidSize=args[x,1], deepSplit=args[x,2])
 	}
 	else{
 		ans<-NULL
 		cat("Nothing fits criteria\n")
 		opt_params<-list(minMidSize=NA, deepSplit=NA)
 	}
 	ans<-ans[order(ans$group),]

 	args<-as.list(match.call())
 	
 	list(sampTab=ans, args=args, opt_params=opt_params)

}


#' butt_dt dynamic tree cluster 
#'
#' butt_dt dynamic tree cluster 
#' @param sampTab sampTab
#' @param datMat data matrix to cluster, cols are samples/cells
#' @param minModSize minimum module size
#' @param deepSplit deepSplit
#' @param distanceMetric distanceMetric
#' @param linkageMetric linkageMetric
#' @param maxTreeHeight maxTreeHeight depends on 
#'
#' @return sampTab with group column indicating cluster assignment
#'
#' @export
#'
.steam_dt<-function
(sampTab,
 hcRes,
 cors,
 minModSize=20,
 deepSplit=FALSE){
 
	###cutResCells<-cutreeDynamicTree(hcRes, minModuleSize=minModSize, deep=deepSplit, maxTreeHeight=maxTreeHeight)
	cutResCells<-cutreeDynamic(hcRes, distM=as.matrix(cors), cutHeight=NULL, minClusterSize=minModSize, deep=deepSplit,  method='hybrid', pamStage=TRUE, pamRespectsDendro = FALSE)
 	cNames<-cutTreeToNames(rownames(sampTab), cutResCells)
	nnames<-names(cNames)
	grp<-rep(0, nrow(sampTab))
	for(nname in nnames){
	  x<-cNames[[nname]]
	  grp[match(x,rownames(sampTab))]<-nname
	}
	sampTabNew<-cbind(sampTab,  group=grp)
	#sampTabNew<-sampTabNew[order(sampTabNew$dtGrp),]

#	list(sampTab=sampTabNew, tree=hcRes, distM=cors)
	sampTabNew
}

#' @export
steam_mclust<-function
(sampTab,
 datMat,
 G=2:9){
	mcRes<-Mclust(datMat, G=G)
	ans<-cbind(sampTab,group=as.character(mcRes$classification))
	ans<-ans[order(ans$group),]
	args<-as.list(match.call())
	opt_params<-list("modelName"=mcRes$modelName, "G"=mcRes$G)
 	list(sampTab=ans, args=args, opt_params=opt_params)
}


# notet aht this will remove group==0 samples
#' @export
steam_dbscan<-function
(sampTab,
 datMat,
 nClusters=c(5,20), # number of clusters
 sClusters=c(10,NA), # size of clusters
 epsRange=seq(0.5,3, by=.1),
 minPtsRange=seq(2,100,by=1),
 seleCrit="median"
 ){
 	if(is.na(sClusters[2])){
 		sClusters[2]<-ceiling(nrow(sampTab)*.75)
 	}
 	cat(sClusters[1],"<->",sClusters[2],"\n")

 	args<-expand.grid(epsRange, minPtsRange)
    nClus<-rep(0, nrow(args))
    minSizes<-rep(0, nrow(args))
    maxSizes<-rep(0, nrow(args))
 	for(i in seq(nrow(args))){
 	##	cat(args[i,1], "::",args[i,2],"\n")
 		tmpAns<-.steam_dbscan(sampTab, datMat, eps=args[i,1],minPts=args[i,2])
 		grps<-table(tmpAns$group)
 		nClus[i]<-length(grps)
 		minSizes[i]<-min(grps)
 		maxSizes[i]<-max(grps)
 	}
 	ans<-data.frame(num_groups=nClus, smallest=minSizes, largest=maxSizes)
 	tmpAns<-cbind(args, ans)
# 	rownames(tmpAns)<-as.character(1:nrow(tmpAns))


 	# select the ones that fit criteria
 	#crit1<-rownames(tmpAns[tmpAns$num_groups> nClusters[1] & tmpAns$num_groups< nClusters[2],])
 	crit1<-which(tmpAns$num_groups> nClusters[1] & tmpAns$num_groups< nClusters[2])
 	crit2<-which(tmpAns$smallest>= sClusters[1] & tmpAns$largest<= sClusters[2])
 	
 	x1<-intersect(crit1, crit2)
 

 	if(length(x1)>=1){


 		if(seleCrit=='median'){
	 		ng_range<-tmpAns[x1,]$num_groups
### 			aa<-median(ng_range)
			aa<-ng_range[floor(length(ng_range)/2)]
 			ai<-which(tmpAns[x1,]$num_groups==aa)[1]
 		}
 		else{
 			ai<-which.min(tmpAns[x1,]$num_groups)
 		}
 		#x<-which.min(tmpAns[x1[ai],]$largets)

 		x<-x1[ai]
 		if(length(x)>1){
 			x<-x[1]
 		}
 		cat("X=", x,"\n")
 		ans<-.steam_dbscan(sampTab, datMat, eps=args[x,1], minPts=args[x,2])
 		opt_params<-list(eps=args[x,1], minPts=args[x,2])
 		cat("eps=",args[x,1], "minPts=",args[x,2],"\n")
 	}
 	else{
 		ans<-NULL
 		cat("Nothing fits criteria\n")
 		opt_params<-list(eps=NA, minPts=NA)
 	}
 	ans<-ans[order(ans$group),]
 	args<-as.list(match.call())
 	
 	ans<-ans[ans$group!=0,]
 	ans<-droplevels(ans)
 	list(sampTab=ans, args=args, opt_params=opt_params)
 
}

#' @export
.steam_dbscan<-function
(sampTab,
 datMat,
 minPts=5,
 eps=1.5){
   dbRes<-dbscan(datMat, eps=eps, minPts=minPts)
   cbind(sampTab, group=as.character(dbRes$cluster))
}





