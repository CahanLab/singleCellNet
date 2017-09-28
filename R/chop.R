# Patrick Cahan (C) 2017
# patrick.cahan@gmail.com

#' @export
gpa<-function
(expDat,
 minClusterSize=20,
 nPCs=2,
 sdExpl=0.01,
 dThresh=0,
 zThresh=2,
 meanType="overall_mean",
 gMax=5){

  # (1) gene stats
  # (2) PCA
  # (3) mclust

  cat("Calculating gene statistics...\n")
  geneStats<-sc_statTab(expDat, dThresh=dThresh)
  cat("PCA...\n")
  pcaRes<-vg_pca(expDat, geneStats, zThresh=zThresh, meanType=meanType)
  ssdds<-pcaRes$pcaRes$sdev
  ans<-0
  if( sum(ssdds[1:nPCs]) / sum(ssdds) > sdExpl){
    cat("MClust...\n")
    res<-simple_steam_mclust(pcaRes$pcaRes$x[,1:nPCs],G=2:gMax)
    sizes<-table(res)
    if(min(sizes)>=minClusterSize){
      ans<-res
    }
  }
  groups<-ans
 # diffExp<-par_findSpecGenes(expDat[pcaRes$varGenes,], db_steamed$steamed$sampTab))
  list(gs=geneStats, pcaRes=pcaRes, groups=groups)
}


#' @export
gpa_recurse<-function
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
    tmpAns<-gpa_break_tree(expAll, grps, isDone=isDone,nPCs=nPCs, minClusterSize=minClusterSize,sdExpl=sdExpl, dThresh=dThresh, zThresh=zThresh, meanType=meanType, nGrps=nGrps, method=method)
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
gpa_break_tree<-function
(expDat,
  grps, # vector of cluster labels
  isDone, # bool vector indicating whether corresponding groups should be subjected to clustering
 ### aveSilhs, # named list of average sihls per group
  minClusterSize=20,
  nPCs=2,
  sdExpl=0.01,
   dThresh=0,
  zThresh=2,
  meanType="overall_mean",
  nGrps=2,
  method="prcomp"
 ){
  ans_grp<-grps
  ans_list<-list()
  uniGrps<-unique(grps)
  dontCluster<-unique(grps[which(isDone==TRUE)])
  uniqGrps<-setdiff(uniGrps, dontCluster)

  pcs<-matrix(0, nrow=ncol(expDat), ncol=2)
  rownames(pcs)<-colnames(expDat)

  for(i in 1:length(uniGrps)){
    uniGrp<-uniGrps[i]
    xi<-which(grps==uniGrp)
    if(length(xi)>=minClusterSize){      
      cat(uniGrp,":" ,length(xi),"\n")
    ###ans_list[[uniGrp]]<-gpa(expDat[,xi], minClusterSize=minClusterSize,nPCs=nPCs,sdExpl=sdExpl,dThresh=dThresh,zThresh=zThresh,meanType=meanType,gMax=gMax)
     
     ##**## ans_list[[uniGrp]]<-gpa_tree(expDat[,xi],nPCs=nPCs,sdExpl=sdExpl,dThresh=dThresh,zThresh=zThresh,meanType=meanType,nGrps=nGrps, method=method)

      gpaRes<-gpa_tree(expDat[,xi],nPCs=nPCs,sdExpl=sdExpl,dThresh=dThresh,zThresh=zThresh,meanType=meanType,nGrps=nGrps, method=method)  
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
  list(groups=ans_grp, res=ans_list, isDone=isDone, pcs=pcs)#, aveSilhs=aveSils)
}


#' @export
gpa_tree<-function
(expDat,
 nPCs=2,
 sdExpl=0.01,
 dThresh=0,
 zThresh=2,
 meanType="overall_mean",
 nGrps=2,
 method="prcomp"){

  # (1) gene stats
  # (2) PCA
  # (3) hclust -> cuttree

  cat("Calculating gene statistics...\n")
  geneStats<-sc_statTab(expDat, dThresh=dThresh)
  cat("PCA...\n")
  pcaRes<-vg_pca(expDat, geneStats, zThresh=zThresh, meanType=meanType, method=method)
  ans<-0
  ans<-simple_steam_cuttree(pcaRes$pcaRes$x[,1:nPCs],nClusters=nGrps)
  groups<-ans$grps
 # diffExp<-par_findSpecGenes(expDat[pcaRes$varGenes,], db_steamed$steamed$sampTab))
  pcs<-pcaRes$pcaRes$x[,1:nPCs]
  rownames(pcs)<-colnames(expDat)
  list(gs=geneStats, pcaRes=pcaRes, groups=groups, pcs=pcs, avesilh=ans$avesilh)
}


#' @export
gpa_db<-function
(expDat,
 minClusterSize=20,
 nPCs=2,
 sdExpl=0.01,
 dThresh=0,
 zThresh=2,
 meanType="overall_mean",
 eps=1,
 minPts=10){

  # (1) gene stats
  # (2) PCA
  # (3) dbscan

  cat("Calculating gene statistics...\n")
  geneStats<-sc_statTab(expDat, dThresh=dThresh)
  cat("PCA...\n")
  pcaRes<-vg_pca(expDat, geneStats, zThresh=zThresh, meanType=meanType)
  ssdds<-pcaRes$pcaRes$sdev
  ans<-0
  if( sum(ssdds[1:nPCs]) / sum(ssdds) > sdExpl){
    cat("DBscan...\n")
    res<-simple_steam_dbscan(pcaRes$pcaRes$x[,1:nPCs],eps=eps, minPts=minPts)
    sizes<-table(res)
    ##if(min(sizes)>=minClusterSize){
      ans<-res
    ##}
  }
  groups<-ans
 # diffExp<-par_findSpecGenes(expDat[pcaRes$varGenes,], db_steamed$steamed$sampTab))
  list(gs=geneStats, pcaRes=pcaRes, groups=groups)
}


#' @export
gpa_break<-function
(expDat,
  grps,
  minClusterSize=20,
  nPCs=2,
  sdExpl=0.01,
 dThresh=0,
 zThresh=2,
 meanType="overall_mean",
 gMax=5){
  ans_grp<-grps
  ans_list<-list()
  uniGrps<-unique(grps)
  for(uniGrp in uniGrps){
    xi<-which(grps==uniGrp)
    ans_list[[uniGrp]]<-gpa(expDat[,xi], minClusterSize=minClusterSize,nPCs=nPCs,sdExpl=sdExpl,dThresh=dThresh,zThresh=zThresh,meanType=meanType,gMax=gMax)
    ans_grp[xi]<-ans_list[[uniGrp]]$groups
  }
  list(groups=ans_grp, res=ans_list)
}


#' @export
gpa_break_db<-function
(expDat,
  grps,
  nPCs=2,
  sdExpl=0.01,
 dThresh=0,
 zThresh=2,
 meanType="overall_mean",
 eps=1,
 minPts=10){
  ans_grp<-grps
  ans_list<-list()
  uniGrps<-unique(grps)
  for(uniGrp in uniGrps){
    xi<-which(grps==uniGrp)
    ###ans_list[[uniGrp]]<-gpa(expDat[,xi], minClusterSize=minClusterSize,nPCs=nPCs,sdExpl=sdExpl,dThresh=dThresh,zThresh=zThresh,meanType=meanType,gMax=gMax)
    ans_list[[uniGrp]]<-gpa_db(expDat[,xi],nPCs=nPCs,sdExpl=sdExpl,dThresh=dThresh,zThresh=zThresh,meanType=meanType,eps=eps,minPts=minPts)
    ans_grp[xi]<-ans_list[[uniGrp]]$groups
  }
  list(groups=ans_grp, res=ans_list)
}




#' @export
chop_tsne<-function
(datMat,
 perplexity,
 theta)
{
  chopped<-to_tsne(datMat, perplexity=perplexity, theta=theta)
  args<-as.list(match.call())
  list(choppedDat=chopped, args=args)
}

#' @export
to_tsne<-function
(datMat, #cols are vars, rows are cells/samples
 perplexity=30,
 theta=0.30){
  tres<-Rtsne(datMat, pca=FALSE, perplexity=perplexity, theta=theta)
  xres<-tres$Y
  colnames(xres)<-c("TSNE.1", "TSNE.2")
  rownames(xres)<-rownames(datMat)
  xres
}

#' @export
chop_pca<-function
(expDat, # rows=vars, cols=cells
 geneStats,
 zThresh,
 meanType,
 varGenes=NULL)
{
  if(is.null(varGenes)){
    vgPCA_res<-vg_pca(expDat, geneStats, zThresh, meanType)
  }
  else{
    vgPCA_res<-list()
    vgPCA_res[['varGenes']]<-varGenes
    vgPCA_res[['pcaRes']]<-prcomp(t(expDat[varGenes,]),center=T,scale=TRUE)
  }
  args<-as.list(match.call())
  list(choppedDat=vgPCA_res[['pcaRes']]$x,  args=args, pcaSD=vgPCA_res[['pcaRes']]$sd, varGenes=vgPCA_res[['varGenes']])
}



#' find var genes and pca
#'
#' find var genes and pca
#'
#' @param expDat expDat
#' @param geneStats, from running sc_statTab
#' @param zThresh thresh for defining genes as variable
#'
#' @return list of varGenes and pcaRes
#' 
#' @export
#'
#' @export
vg_pca<-function
(expDat,
 geneStats,
 zThresh=2,
 meanType="overall_mean",
 method='prcomp',
 max.iter=10)
{
  ans<-list()
  ans[['varGenes']]<-findVarGenes(expDat,geneStats,zThresh=zThresh, meanType=meanType)
  if(method=='prcomp'){
    ans[['pcaRes']]<-prcomp(t(expDat[ ans[['varGenes']],]),center=T,scale=TRUE)
  }
  else{
    # will not have a sdev item
    normDat<-prep(t(expDat[ans[['varGenes']],]), scale="uv", center=TRUE)
    tmpAns<-rpca(normDat, trace=FALSE, max.iter=max.iter)
    tmpX<-tmpAns$L.svd$u
    rownames(tmpX)<-colnames(expDat)
    ans[['pcaRes']]<-list(x=tmpX)
  }
  ans
}

#' find topX PCs
#'
#' find topX PCs
#'
#' @param vp result of running vg_pca
#' @param qtile quantile
#'
#' @return index of PC at which _quantile_ of sd is included
#' 
#' @export
find_top_pc<-function
(vp,
 qtile=.25)
{
  x<-quantile(vp$pcaRes$sd, qtile)
  max(which(vp$pcaRes$sd>x))
}

#' run dbscan, remove bad points, run again, plot 2 plots, return df of dbscan res
#'
#' run dbscan, remove bad points, run again, plot 2 plots, return df of dbscan res
#'
#' @param mtsneRes, # has columns TSNE.1 and TSNE.1
#' @param eps dbscan eps
#' @param minPts dbscan minPts 
#'
#' @return ggplot
#' 
#' @export
#'
pruneAndScan<-function 
(tsneRes, # has columns TSNE.1 and TSNE.1
 eps=1.8,
 minPts=10){
  ans<-list()
  # remove any linger 'group' columns
  cnames<-colnames(tsneRes)
  cnames<-setdiff(cnames, "group")
  tsneRes<-tsneRes[,cnames]

  dbRes<-dbscan(tsneRes[,c("TSNE.1", "TSNE.2")], eps=eps, minPts=minPts)
  tmpAns<-cbind(tsneRes, group=as.character(dbRes$cluster))
  plot1<-plotDBscan(tmpAns)
  prunedRes<-tmpAns[which(tmpAns$group!=0),]
  plot2<-plotDBscan(prunedRes)
  multiplot(plot1, plot2, cols=2)
  tmpAns
}
