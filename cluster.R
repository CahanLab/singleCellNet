

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
vg_pca<-function
(expDat,
 geneStats,
 zThresh=2,
 meanType="overall_mean")
{
  ans<-list()
  ans[['varGenes']]<-findVarGenes(expDat,geneStats,zThresh=zThresh, meanType=meanType)
  ans[['pcaRes']]<-prcomp(t(expDat[ ans[['varGenes']],]),center=T,scale=TRUE)
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
#'
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

cutTreeToGenes<-function # converts the result of running cutreeDynamicTree to a list of genes
(genes, # vector of gene names in same order as ctRes
 ctRes # result of running cutreeDynamicTree 
){
  ans<-list()
  grpNames<-unique(ctRes)
  for(grpName in grpNames){
    ans[[as.character(grpName)]]<-sort(genes[ which(ctRes==grpName) ])
  }
  ans;
}