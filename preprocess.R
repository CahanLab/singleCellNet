# scCellNet
# (C) Patrick Cahan 2012-2017

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

#' weighted subtraction from mapped reades and log applied to all
#'
#' Simulate expression profile of  _total_ mapped reads
#' @param expRaw matrix of total mapped reads per gene/transcript
#' @param total numeric post transformation sum of read counts
#'
#' @return vector of downsampled read mapped to genes/transcripts
#'
#' @export
trans_rnaseq<-function
(expRaw,
 total
 ){
    expCountDnW<-apply(expRaw, 2, downSampleW, total)
    log(1+expCountDnW)
  }

  #' express as a fraction of the column total
#'
#' express as a fraction of the column total
#' @param expDat expression matrix
#' @param xFact scale by this value
#'
#' @return transformed data matrix
trans_dNorm<-function 
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
  ans;
}

#' shift the expression of each sample so that the min=minVal
#'
#' shift the expression of each sample so that the min=minVal
#' @param expDat expression matrix
#' @param minVal
#'
#' @return transformed expression matrix
trans_eShift<-function #
(expDat,
 minVal=0){
  
  myFunc<-function(vect, myMin){
    vect-min(vect)+myMin;
  }
  
  ans<-apply(expDat, 2, myFunc, myMin=minVal);
  rownames(ans)<-rownames(expDat);
  colnames(ans)<-colnames(expDat);
  ans;
}

#' quantile normalize expression matrix
#'
#' equalize the distributions of columns
#' @param expDat
#'
#' @return quantile normalized matrix
Norm_quantNorm<-function#
(expDat){
  require(preprocessCore);
  if(is.data.frame(expDat)){
    expDat<-as.matrix(expDat);
  }
  ans<-normalize.quantiles(expDat);
  colnames(ans)<-colnames(expDat);
  rownames(ans)<-rownames(expDat);
  ans;
}