# Patrick Cahan (C) 2017
# patrick.cahan@gmail.com

#'export




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
