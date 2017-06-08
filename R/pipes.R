# singleCellNet
# (C) Patrick Cahan 2012-2017



pipe_averageGroup<-function
(washed,
pipeSteamed
){

	varGenes<-pipeSteamed[['cp_pca']][['varGenes']]
	stx<-pipeSteamed[['steamed']][['sampTab']]
	grps<-unique(stx$group)
	expDat<-washed[['expDat']]
	ans<-matrix(0, nrow=length(varGenes), ncol=length(grps))
	colnames(ans)<-grps
	rownames(ans)<-varGenes
	for(grp in grps){
		stTmp<-stx[stx$group==grp,]
		expTmp<-expDat[varGenes,rownames(stTmp)]
		ans[,grp]<-apply(expTmp,1, mean)
	}
	ans
}


# splits data
# makes classifiers
# apply to held out data
# returns classifier too
#' @export
pipe_butter<-function
(pipeSteamed, # has cpName(choppedDat), varGenes
 washedDat, # has expDat
 propTrain=0.25,
 partOn="group",
 nTrees=200)
{
  geneLists<-list()
  stDat<-pipeSteamed[['steamed']][['sampTab']]
  predictors<-pipeSteamed[['cp_pca']][['varGenes']]
  cts<-as.vector(unique(stDat[,partOn]))
  for(ct in cts){
    geneLists[[ct]]<-predictors
  }

  # split into training and test data
  ttList<-divide_sampTab(stDat, propTrain, partOn)
  stTrain<-ttList[['stTrain']]

  # make RFs
  myRFs<-makeRFs(washedDat[['expDat']][predictors,rownames(stTrain)], stTrain, geneLists, dLevel=partOn,
    nTrees=nTrees)
  
  # classify held out data
  stVal<-ttList[['stVal']]

  list(classResVal=list(classRes=sc_classify(myRFs, washedDat[['expDat']][predictors,rownames(stVal)], geneLists),
  						stVal=stVal),
  		classifiers=myRFs,
  		predictors=geneLists)
}

# chop and steam dbscan style
#' @export
pipe_dbscan<-function
(washedDat,
 sampTab,
 topPC,
 nClusters=c(5,20),
 sClusters=c(10,NA),
 zThresh=2)
{
	cp_pca<-chop_pca(washedDat[['expDat']], washedDat[['geneStats']], 
		zThresh=zThresh, meanType="overall_mean")

	# tsne
	cat("tsne-ing\n")
	cp_tsne<-chop_tsne(cp_pca[['choppedDat']][,1:topPC], perplexity=30, theta=.3)

	# steam  dbscan	
	cat("dbscan\n")
	stm_dbs<-steam_dbscan(sampTab, cp_tsne[['choppedDat']][,1:2],nClusters=nClusters,sClusters=sClusters)

	list(cp_pca=cp_pca, cp_tsne=cp_tsne, steamed=stm_dbs)
}

# chop and steam dynamiccuttree style
#' @export
pipe_treecut<-function
(washedDat,
 sampTab,
 topPC,
 nClusters=c(5,20),
 sClusters=c(10,NA),
 zThresh=2)
{
	cp_pca<-chop_pca(washedDat[['expDat']], washedDat[['geneStats']], zThresh=zThresh, meanType="overall_mean")

	# steam  dbscan	
	cat("cuttree\n")
	stm<-steam_treecut(sampTab, cp_pca[['choppedDat']][,1:topPC],nClusters=nClusters,sClusters=sClusters)

	list(cp_pca=cp_pca, steamed=stm)
}


# input: chopped data
# out: df of cAss with steam_method column
#' @export
pipe_cAss<-function
(washedDat, # expDat, geneStats
 sampTab, 
 topPC=20){

	methods<-c("dbscan", "dynamicTreeCut", "mclust")
	mRes<-list()

	# all methods need PCA dimension reduction, this will gives us varGenes
	cat("reducing dimensionality\n")
	cp_pca<-chop_pca(washedDat[['expDat']], washedDat[['geneStats']], 
		zThresh=2, meanType="overall_mean")

	# dbscan needs tsne
	cat("tsne-ing\n")
	cp_tsne<-chop_tsne(cp_pca[['choppedDat']][,1:topPC], perplexity=30, theta=.3)

	cat("Mclust\n")
	# steam MClust
	stm_mc<-steam_mclust(sampTab,cp_pca[['choppedDat']][,1:3], G=4:15)

	# steam  dbscan	
	cat("dbscan\n")
	stm_dbs<-steam_dbscan(sampTab, cp_tsne[['choppedDat']][,1:2])

	# steam treecut
	cat("dtree cut\n")
	stm_tree<-steam_treecut(sampTab, cp_pca[['choppedDat']][,1:topPC])


	# pre-butter each
	cat("pre butter manual\n")	
	classRes_manual<-prebutter(cp_pca$varGenes, washedDat[['expDat']], sampTab,
		propTrain=.25, partOn="prefix", nTrees=200)

	cat("pre butter dbscan\n")
	classRes_db<-prebutter(cp_pca$varGenes, washedDat[['expDat']], stm_dbs[['sampTab']],
		propTrain=.25, partOn="group", nTrees=200)
	cat("pre butter mclust\n")	
	classRes_mclust<-prebutter(cp_pca$varGenes, washedDat[['expDat']], stm_mc[['sampTab']],
		propTrain=.25, partOn="group", nTrees=200)
	cat("pre butter tree\n")	
	classRes_treeCut<-prebutter(cp_pca$varGenes, washedDat[['expDat']], stm_tree[['sampTab']],
		propTrain=.25, partOn="group", nTrees=200)
	

	# assess
	cAss_manual<-easy_assess(classRes_manual[['classRes']], classRes_manual[['stVal']])
	cAss_manual<-cbind(cAss_manual, group=cAss_manual[,'prefix'])
	cAss_dbscan<-easy_assess(classRes_db[['classRes']], classRes_db[['stVal']])
	cAss_mclust<-easy_assess(classRes_mclust[['classRes']], classRes_mclust[['stVal']])
	cAss_treeCut<-easy_assess(classRes_treeCut[['classRes']], classRes_treeCut[['stVal']])

	ans<-cbind(cAss_dbscan, method=rep("dbscan", nrow(cAss_dbscan)))
	ans<-rbind(ans, cbind(cAss_mclust, method=rep("mclust", nrow(cAss_mclust))))
	ans<-rbind(ans, cbind(cAss_treeCut, method=rep("dtreeCut", nrow(cAss_treeCut))))
	ans<-rbind(ans, cbind(cAss_manual, method=rep("manual", nrow(cAss_manual))))
	ans
}

#' @export
pipe_steam_list<-function
(washedDat,
 sampTab,
 topPC=20){

	# pca
	cp_pca<-chop_pca(washedDat[['expDat']], washedDat[['geneStats']], 
		zThresh=2, meanType="overall_mean")

	# dbscan needs tsne
	cat("tsne-ing\n")
	cp_tsne<-chop_tsne(cp_pca[['choppedDat']][,1:topPC], perplexity=30, theta=.3)

	# steam  dbscan	
	cat("dbscan\n")
	stm_dbs<-steam_dbscan(sampTab, cp_tsne[['choppedDat']][,1:2])

	cat("Mclust\n")
	# steam MClust
###	stm_mc<-steam_mclust(sampTab,cp_pca[['choppedDat']][,1:3], G=4:15)
	stm_mc<-steam_mclust(sampTab,cp_tsne[['choppedDat']][,1:2], G=4:15)

	# steam treecut
	cat("dtree cut\n")
	stm_tree<-steam_treecut(sampTab, cp_pca[['choppedDat']][,1:topPC])

	list(cp_pca=cp_pca, cp_tsne=cp_tsne, 
		manual=sampTab, 
		dbscan=stm_dbs[['sampTab']],
		mc=stm_mc[['sampTab']],
		dtree=stm_tree[['sampTab']])
}

# input:
# chopped data 
# expDat (transformed as you wish)
# sampTab
# steam method
# steamParams (list)
# out: df of cAss with steam_method column
#' @export
pipe_cAss_all<-function
(steamed, # result of running `pipe_steam_list`
 expDat, # expression data as you like it to be transformed
 sampTab){
	
	cp_pca<-steamed[['cp_pca']]
	cp_tsne<-steamed[['cp_tsne']]

	# pre-butter each
	cat("pre butter manual\n")	
	classRes_manual<-prebutter(cp_pca$varGenes, expDat, sampTab,
		propTrain=.25, partOn="prefix", nTrees=200)

	cat("pre butter dbscan\n")
	classRes_db<-prebutter(cp_pca$varGenes, expDat, steamed[['dbscan']],
		propTrain=.25, partOn="group", nTrees=200)

	cat("pre butter mclust\n")	
	classRes_mclust<-prebutter(cp_pca$varGenes, expDat, steamed[['mc']],
		propTrain=.25, partOn="group", nTrees=200)

	cat("pre butter tree\n")	
	classRes_treeCut<-prebutter(cp_pca$varGenes, expDat, steamed[['dtree']],
		propTrain=.25, partOn="group", nTrees=200)

	# assess# assess
	cAss_manual<-easy_assess(classRes_manual[['classRes']], classRes_manual[['stVal']])
	cAss_manual<-cbind(cAss_manual, group=cAss_manual[,'prefix'])
	cAss_dbscan<-easy_assess(classRes_db[['classRes']], classRes_db[['stVal']])
	cAss_mclust<-easy_assess(classRes_mclust[['classRes']], classRes_mclust[['stVal']])
	cAss_treeCut<-easy_assess(classRes_treeCut[['classRes']], classRes_treeCut[['stVal']])

	ans<-cbind(cAss_dbscan, method=rep("dbscan", nrow(cAss_dbscan)))
	ans<-rbind(ans, cbind(cAss_mclust, method=rep("mclust", nrow(cAss_mclust))))
	ans<-rbind(ans, cbind(cAss_treeCut, method=rep("dtreeCut", nrow(cAss_treeCut))))
	ans<-rbind(ans, cbind(cAss_manual, method=rep("manual", nrow(cAss_manual))))
	ans
}






