# singleCellNet
# (C) Patrick Cahan 2012-2017


#' subsets on ok genes, applies selected data transform, and recalculates gene stats
#'
#' down samples, computes per gene stats, finds ok genes
#' @param preWashed list of expDat=expDat, geneStats=geneStats, okGenes=okGenes from prewash()
#' @param transMethod transMethod 
#' @param threshold threshold if method==binary then values less than threshold get set to 0, everything else gets set to 1
#' @param removeBad
#'
#' @return list of geneStats, expDat, transMethod
#'
#' @export
wash<-function
(preWashed,
 transMethod="prop",
 threshold=1,
 removeBad=TRUE){

    geneStats<-preWashed[['geneStats']]
    expDat<-preWashed[['expDat']]

    if(removeBad){     
      geneStats<-geneStats[preWashed[['okGenes']],]
      expDat<-expDat[preWashed[['okGenes']],]
    }

    if(transMethod=="binary"){
      expDat<-trans_binarize(expDat, threshold=1)
    }
    else{
      if(transMethod=="zscore"){
        expDat<-trans_zscore(expDat)
      }
      else{
        expDat<-trans_prop(expDat)
      }
    }

    geneStats<-sc_statTab(expDat, dThresh=0)
    list(geneStats=geneStats, expDat=expDat, transMethod=transMethod)
}

#' down samples, computes per gene stats, finds ok genes
#'
#' down samples, computes per gene stats, finds ok genes
#' @param expDat expDat
#' @param stDat sample table
#' @param countDepth counts to sample per cell
#' @param minGenes min genes 
#' @param alpha1 alpha 1
#' @param mu mu 
#' @param thresh detection threshold for gene stats
#'
#' @return list of expDat, geneStats, okgenes
#'
#' @export
prewash<-function
(expDat,
 stDat,
 countDepth=1e3,
 minGenes=10,
 alpha1=0.01,
 mu=1,
 thresh=0){

  expDat<-weighted_down(expDat,countDepth)
  geneStats<-sc_statTab(expDat, dThresh=thresh)
  alpha2<-minGenes/nrow(stDat)
  okGenes<-sc_filterGenes(geneStats, alpha1=alpha1, alpha2=alpha2, mu=mu)
  list(expDat=expDat, geneStats=geneStats, okGenes=okGenes)
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
 total=1e5,
 dThresh=0){ 

  totalSignal<-sum(vector)
  wAve<-vector/totalSignal
###  resid<-sum(vector)-total #num to subtract from sample
  resid<-totalSignal-total #num to subtract from sample
  residW<-wAve*resid # amount to substract from each gene
  ans<-vector-residW
  ans[which(ans<dThresh)]<-0
  ans
}



# not used
newWD<-function(expDat, 
  total=1e5, 
  dThresh=0){
  cSums  <- colSums(expDat)
  resids <- cSums - total
  props  <- sweep(expDat, MARGIN=2, FUN="/", STATS=colSums(expDat))
  
  rrids  <- cSums - total

  # tmpAns <- expDat - t(t(props) * rrids)
  tmpAns <- expDat - sweep(props, 2, rrids, "*")
  tmpAns[which(tmpAns<dThresh)] <- 0
  tmpAns

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
(expDat,
 total,
 dThresh=0
 ){
  if(class(expDat)[1]!='matrix'){
    cSums  <- Matrix::colSums(expDat)
    props <- Matrix::t(expDat) / cSums
    rrids  <- cSums - total
    tmpAns <- expDat - Matrix::t(props * rrids)
    tmpAns[Matrix::which(tmpAns<dThresh)] <- 0
  }
  else{
    cSums  <- colSums(expDat)
    props <- t(expDat) / cSums
    rrids  <- cSums - total
    tmpAns <- expDat - t(props * rrids)
    tmpAns[which(tmpAns<dThresh)] <- 0
  }
  
  tmpAns
}




if(FALSE){

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
 total,
 dThresh=0
 ){
    expCountDnW<-apply(expRaw, 2, downSampleW, total=total, dThresh=dThresh)
    #log(1+expCountDnW)
    expCountDnW
  }
}

if(FALSE){

#' @export
trans_prop<-function 
(expDat,
 xFact=1e5
){
  ans<-matrix(0, nrow=nrow(expDat), ncol=ncol(expDat))
  for(i in seq(ncol(expDat))){
    ans[,i]<-expDat[,i]/sum(expDat[,i]);    
  }
  ans<-ans*xFact;
  colnames(ans)<-colnames(expDat);
  rownames(ans)<-rownames(expDat);
  log(1+ans)
}
}

#' @export
trans_prop<-function
(expDat,
  xFact=1e5){

  if(class(expDat)[1]!='matrix'){
    cSums  <- Matrix::colSums(expDat)
    ans <- Matrix::t(log(1 + xFact * Matrix::t(expDat) / cSums))
  }
  else{
    cSums  <- colSums(expDat)
    ans<-  t(log(1 + xFact * t(expDat) / cSums))
  }
  ans
 
}

#' @export
trans_zscore<-function
(expRaw){
  apply(expRaw, 2, scale)
}

#' @export
trans_binarize<-function
(expRaw,
  threshold=1){
  expRaw[expRaw<threshold]<-0
  expRaw[expRaw>0]<-1
  expRaw
}

