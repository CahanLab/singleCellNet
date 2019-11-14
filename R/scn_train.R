#' @title
#' Training
#' @description
#' Tranining broad class classifier
#' @param stTrain a dataframe that matches the samples with category
#' @param expTrain the expression matrix
#' @param dLevel the name of the column that contains categories
#' @param colName_samp the name of the column that contains sample names
#' @param nTopGenes the number of classification genes per category
#' @param nTopGenePairs the number of top gene pairs per category
#' @param nRand number of random profiles generate for training
#' @param nTrees number of trees for random forest classifier
#' @param weightDown_total numeric post transformation sum of read counts for weighted_down function
#' @param weightedDown_dThresh the threshold at which anything lower than that is 0 for weighted_down function
#' @param transprop_xFact scaling factor for transprop
#'
#' @return a list containing normalized expression data, classification gene list, cnPRoc
#' @export

scn_train <- function(stTrain, 
		      expTrain, 
		      dLevel, 
          	      colName_samp="row.names", 
	  	      nTopGenes = 10, 
		      nTopGenePairs = 25, 
		      nRand = 70, 
		      nTrees = 1000,
          	      stratify=FALSE, 
		      weightedDown_total = 1e4, 
 	 	      weightedDown_dThresh = 0.25) {

   if (class(stTrain) != "data.frame") {
      stTrain<-as.data.frame(stTrain)
   }

   if (colName_samp != "row.names") {
     rownames(stTrain)<-stTrain[, colName_samp]
   }

   cat("Sample table has been prepared\n")

   expTnorm<-trans_prop(expTrain, weightedDown_total, dThresh = weightedDown_dThresh)
   cat("Expression data has been normalized\n")

   cat("Finding classification genes\n")
   system.time(cgenes<-findClassyGenes(expDat = expTnorm, sampTab = stTrain, dLevel = dLevel, topX = nTopGenes))

   cgenesA<-cgenes[['cgenes']]
   cat("There are ", length(cgenesA), " classification genes\n")

   grps<-cgenes[['grps']]

   #catch errors when there is NA or emtpy string in cluster/cell type label
   if(sum(grps == "")>1 | sum(is.na(grps))>1){
    stop("There is NA or empty string in your dLevel. Please remove them before proceeding.")
   }
   
   cat("Finding top pairs\n")

   system.time(xpairs<-ptGetTop(expDat = expTrain[cgenesA,], cell_labels = grps, cgenes_list = cgenes[['cgenes_list']], topX=nTopGenePairs, sliceSize=5000))
   cat("There are", length(xpairs), "top gene pairs\n")

   system.time(pdTrain<-query_transform(expTrain[cgenesA, ], xpairs))
   cat("Finished pair transforming the data\n")

   tspRF<-sc_makeClassifier(pdTrain[xpairs,], genes=xpairs, groups=grps, nRand = nRand, ntrees = nTrees, stratify=stratify)
   cnProc<-list("cgenes"= cgenesA, "xpairs"=xpairs, "grps"= grps, "classifier" = tspRF)

   returnList<-list("sampTab" = stTrain, "cgenes_list" = cgenes[['cgenes_list']], "cnProc" = cnProc)

   cat("All Done\n")
   #return
   returnList
}


 
#' Predict query using broad class classifier 
#' @description
#' The function predicts the query data using the broad class classifier 
#' @param cnProc the cnProc of the broad classifier 
#' @param expDat the expression data of query data 
#' @param nrand the number of random profiles generate for evaluation process 
#' @return a classification score matrix 
#'
#' @export 
scn_predict<-function(cnProc, expDat, nrand = 2) {

   expVal = expDat 

   rf_tsp<-cnProc[['classifier']]
   cgenes<-cnProc[['cgenes']]
   xpairs<-cnProc[['xpairs']]

   cat("Loaded in the cnProc\n")

   expValTrans<-query_transform(expVal[cgenes,], xpairs)
   classRes_val<-rf_classPredict(rf_tsp, expValTrans, numRand=nrand)

   cat("All Done\n")

   #return 
   classRes_val
}

#' transform mclust class res to expression prob estimates
#'
#' transform mclust class res to expression prob estimates
#' @param mcRes what is returned by running Mclust
#' @param expDat matrix used as input to Mclust
#' @param pres bool whehter to report prob of expression highest 2 classes or not only highest class
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
#' @export
#'
mclust_Mat<-function
(expDat,
 pres=TRUE){

        cat("starting mclust...\n")
        acounts<-Mclust(as.vector(expDat), G=3)
        cat("done with mclust\n")


        probTrans(acounts, expDat, pres=pres)
}

#' harmonize data sets
#' @param expList list of expression matricies
#' @param dsVal target expression count
#' @param pres bool whehter to report prob of expression highest 2 classes; default or not only highest class
#'
#' @return list1=list of normalized expression matricies, list2=list of prob of exp matricies
#' @export
#'
harmonize<-function(
 expList,
 trans,#vector of booleans indicating whether to downsample prior to prob
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

#' @title
#' Make Classifier
#'
#' @description
#' Create a random forest classifier with the transformed training data from \code{\link{query_transform}}.
#'
#' @param expTrain transformed training data from \code{\link{query_transform}}
#' @param genes vector of gene pairs from \code{\link{ptGetTop}} used as predictors
#' @param groups named vector of cells to cancer categories
#' @param nRand number of randomized profiles to make
#' @param ntrees number of trees to build
#' @param stratify whether to use stratified sampling or not
#' @param samplesize the samplesize for straified sampling
#' @importFrom randomForest randomForest
#'
#' @return Random Forest Classifier object
#' @export
sc_makeClassifier<-function(expTrain, genes, groups, nRand=70, ntrees=2000, stratify=FALSE, sampsize=40){
  randDat<-randomize(expTrain, num=nRand)
  #randDat<-ModifiedRandomize(expTrain, num=nRand)

  expTrain<-cbind(expTrain, randDat)

  allgenes<-rownames(expTrain)
  missingGenes<-setdiff(unique(genes), allgenes)
  cat("Number of missing genes ", length(missingGenes),"\n")
  ggenes<-intersect(unique(genes), allgenes)

  # return random forest object
  if(!stratify){
    randomForest::randomForest(t(expTrain[ggenes,]), as.factor(c(groups, rep("rand", ncol(randDat)))), ntree=ntrees)
  }else{
    randomForest::randomForest(t(expTrain[ggenes,]), as.factor(c(groups, rep("rand", ncol(randDat)))), ntree=ntrees, strata = as.factor(c(groups, rep("rand", ncol(randDat)))), sampsize=rep(sampsize, length(c(unique(groups), "rand"))))
  }
}

#' classify samples
#'
#' classify samples
#' @param rfObj result of running sc_makeClassifier
#' @param expQuery expQuery
#' @param numRand numRand

#' @return classRes matrix
#' @export
#'
rf_classPredict<-function(
  rfObj,
  expQuery,
  numRand=50){
      if(numRand > 0 ){
        randDat<-randomize(expQuery, num=numRand)
        expQuery<-cbind(expQuery, randDat)
      }
   	preds<-rownames(rfObj$importance)
        xpreds<-t(predict(rfObj, t(expQuery[preds,]), type='prob'))
        colnames(xpreds)<-colnames(expQuery)
        xpreds
}

