# scCellNet
# (C) Patrick Cahan 2012-2017

wash<-function
(expDat,
  prewashed){
    geneStats<-pwashed[['geneStats']][pwashed[['okGenes']],]
    expDat<-expDat[pwashed[['okGenes']],]
    list(geneStats=geneStats,expDat=expDat)
}

prewash<-function
(expDat,
 stDat,
 minGenes=10,
 alpha1=0.01,
 mu=1,
 thresh=0){

  geneStats<-sc_statTab(expDat, dThresh=thresh)
  alpha2<-minGenes/nrow(stDat)
  okGenes<-sc_filterGenes(geneStats, alpha1=alpha1, alpha2=alpha2, mu=mu)
  list(geneStats=geneStats, okGenes=okGenes)
}



#' weighted subtraction from mapped reades
#'
#' Simulate expression profile of  _total_ mapped reads
#' @param vector of total mapped reads per gene/transcript
#' @param total post transformation sum of read counts
#'
#' @return vector of downsampled read mapped to genes/transcripts
#'
#' @export
downSampleW<-function
(vector,
total=1e5){ 

  totalSignal<-sum(vector)
  wAve<-vector/totalSignal
  resid<-sum(vector)-total #num to subtract from sample
  residW<-wAve*resid # amount to substract from each gene
  ans<-vector-residW
  ans[which(ans<0)]<-0
  ans
}

#' weighted subtraction from mapped reades, applied to all
#'
#' Simulate expression profile of  _total_ mapped reads
#' @param expRaw matrix of total mapped reads per gene/transcript
#' @param total numeric post transformation sum of read counts
#'
#' @return vector of downsampled read mapped to genes/transcripts
#'
#' @export
weighted_down<-function
(expRaw,
 total
 ){
    expCountDnW<-apply(expRaw, 2, downSampleW, total)
    #log(1+expCountDnW)
    expCountDnW
  }

trans_prop<-function 
(expDat,
 xFact=1e5
){
  ans<-matrix(0, nrow=nrow(expDat), ncol=ncol(expDat));
  for(i in seq(ncol(expDat))){
    ans[,i]<-expDat[,i]/sum(expDat[,i]);    
  }
  ans<-ans*xFact;
  colnames(ans)<-colnames(expDat);
  rownames(ans)<-rownames(expDat);
  log(1+ans)
}

trans_zscore<-function
(expRaw){
  apply(expRaw, 2, scale)
}

trans_binarize<-function
(expRaw,
  threshold=1){
  expRaw[expRaw<threshold]<-0
  expRaw[expRaw>0]<-1
  expRaw
}

 


  #' express as a fraction of the column total
#'
#' express as a fraction of the column total
#' @param expDat expression matrix
#' @param xFact scale by this value
#'
#' @return transformed data matrix