# Patrick Cahan (C) 2017
# patrick.cahan@gmail.com

#
#' transform mclust class res to expression prob estimates
#'
#' transform mclust class res to expression prob estimates
#' @param mcRes what is returned by running Mclust
#' @param expDat matrix used as input to Mclust
#' @param pres bool whehter to report prob of expression (highest 2 classes) or not (only highest class)
#'
#' @return matrix of probabilities of expression class for each cell and gene
#' @export
#'
probTrans<-function(
  mcRes,
  expDat,
  pres=TRUE
  )
  {

  	# need to determine whether labels correlate with mean expression
  	meansClasses<-rep(0,3)
  	for(i in 1:3){
	  meansClasses[i]<-mean(mcRes$data[which(mcRes$classification==i)])
	}
	cLow<-which.min(meansClasses)
	cHigh<-which.max(meansClasses)


    if(pres){
      za<-1-mcRes$z[,cLow] # prob of detected
    }
    else{
      za<-mcRes$z[,cHigh] #prob high
    }
    zz<-matrix(za, nrow=nrow(expDat))
     colnames(zz)<-colnames(expDat)
     rownames(zz)<-rownames(expDat)
    zz
  }

#' weighted subtraction from mapped reades and log applied to all
#'
#' down samples
#' @param expRaw matrix of total mapped reads per gene/transcript
#' @param total numeric post transformation sum of read counts
#' @param genes subset of genes to perform this on
#'
#' @return vector of downsampled read mapped to genes/transcripts
#'
#' @export
sc_trans_rnaseq<-function
(expRaw,
 total,
 genes
 ){
    expCountDnW<-apply(expRaw[genes,], 2, downSampleW, total)
    log(1+expCountDnW)
  }

#
#' run Mclust on an expression matrix
#'
#' run Mclust on an expression matrix
#' @param expDat matrix
#' @param pres boolean 
#' 
#' @return matrix rows=genes, cols=cells values 
#' @export bool whehter to report prob of expression (highest 2 classes) or not (only highest class)
#'
mclust_Mat<-function
(expDat,
 pres=TRUE){

	cat("starting mclust...\n")
	acounts<-Mclust(as.vector(expDat), G=3)
	cat("done with mclust\n")


	probTrans(acounts, expDat, pres=pres)
}


#
#' harmonize data sets
#'
#' harmonize data sets
#' @param expList list of expression matricies
#' @param dsVal target expression count
#' @param pres bool whehter to report prob of expression (highest 2 classes; default) or not (only highest class)
#'
#' @return list1=list of normalized expression matricies, list2=list of prob of exp matricies
#' @export
#'
harmonize<-function(
 expList,
 dsVal=1e6,
 pres=TRUE){

	# first find the genes in common to all matricies
	allgenes<-lapply(expList, rownames)
	cgenes<-Reduce(intersect, allgenes)

	# transform each matrix
	transList<-lapply(expList, sc_trans_rnaseq, total=dsVal, genes=cgenes)

	# convert to probabilities
	probList<-lapply(transList, mclust_Mat, pres=pres)

	ans<-list(transList=transList, probList=probList)
	class(ans)<-"Harmonized"
	ans
}



#
#' randomize data matrix 
#'
#' randomize data matrix 
#' @param expDat expDat
#' @param num number of profiles to return
#'
#' @return exp matrix random
#' @export
#'
randomize<-function(
 expDat,
 num=50){


	randDat<-t(apply(expDat, 1, sample))	
 	randDat<-apply(randDat, 2, sample)

 	randDat<-randDat[,sample(1:ncol(randDat), num)]
	colnames(randDat)<-paste0(rep("rand_", num), 1:num)
	rownames(randDat)<-rownames(expDat)
	randDat
}

#
#' make a classifier, with a randomized class too
#'
#' rmake a classifier, with a randomized class too
#' @param expTrain training data
#' @param genes vector of genes to use as predictors
#' @param groups named vector of cells to groups or classes
#' @param nRand =50 num of randomized profiles to make
#' @ntrees ntrees =2000 number of trees to build

#' @return exp matrix random
#' @export
#'
sc_makeClassifier<-function(
  expTrain,
  genes,
  groups,
  nRand=50,
  ntrees=2000){


	randDat<-randomize(expTrain, num=nRand)
	expTrain<-cbind(expTrain, randDat)

	allgenes<-rownames(expTrain)

	missingGenes<-setdiff(unique(genes), allgenes)
	cat("Number of mussing genes ", length(missingGenes),"\n")
	ggenes<-intersect(unique(genes), allgenes)
	randomForest(t(expTrain[ggenes,]), as.factor(c(groups, rep("rand", ncol(randDat)))), ntree=2000)

}

rf_classPredict<-function(
  rfObj,
  expQuery,
  numRand=50){

  	randDat<-randomize(expQuery, num=numRand)
  	expQuery<-cbind(expQuery, randDat)

    preds<-rownames(rfObj$importance)
  	xpreds<-t(predict(rfObj, t(expQuery[preds,]), type='prob'))
	colnames(xpreds)<-colnames(expQuery)
	xpreds
}


