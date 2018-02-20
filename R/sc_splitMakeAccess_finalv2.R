
sc_splitMakeAssess<-function
(stDat,
 washedDat,
 xTree2,
 dLevelSID="sample_name",
 propTrain=0.25,
 numGenes = 500,
 nTrees=200){
  
  classifiers <- list()
  geneLists <- list()
  
  topL <- xTree2$grp_list[[2]]
  stDat_topL <- cbind(stDat, topL)
  diffExp <- gnrAll(washedDat[['expDat']], stDat_topL$topL)
  cts<-as.vector(unique(stDat_topL[,"topL"]))
  
  for(i in 1:length(cts)){
    predictors <- getClassGenes(diffExp[[i]],numGenes)
    geneLists[[cts[i]]] <- predictors
  }
  
  
  # split into training and test data
  ttList<-divide_sampTab(stDat_topL, propTrain, dLevel = "topL")
  stTrain<-ttList[['stTrain']]
  
  cat("Making topLevel classifiers\n")
  testRFs<-makeRFs(washedDat[['expDat']][,rownames(stTrain)], stTrain, geneLists, dLevel="topL",
                   nTrees=nTrees)
  
  # classify held out data
  stVal<-ttList[['stVal']]
  ansVal=sc_classify(testRFs, washedDat[['expDat']][,rownames(stVal)], geneLists)
  
  cat("Assessing topLevel...\n")
  assessed <-sc_classAssess(ansVal, stVal, classLevels="topL", dLevelSID=dLevelSID, resolution=0.01)
  classifiers[['topLevel']] <- list(classifiers=testRFs, classRes=ansVal, PRs=assessed, stVal = stVal)
  
  #bottom level
  if(length(xTree2$grp_list) > 2){
    bottomL <- xTree2$grp_list[[length(xTree2$grp_list)]]
    stDat_bottomL <- cbind(stDat, bottomL)
    diffExp <- gnrAll(washedDat[['expDat']], stDat_bottomL$bottomL)
    cts<-as.vector(unique(stDat_bottomL[,"bottomL"]))
    
    for(i in 1:length(cts)){
      predictors <- getClassGenes(diffExp[[i]],numGenes)
      geneLists[[cts[i]]] <- predictors
    }
    # split into training and test data
    ttList<-divide_sampTab(stDat_bottomL, propTrain, dLevel = "bottomL")
    stTrain<-ttList[['stTrain']]
    
    cat("Making bottomLevel classifiers\n")
    testRFs<-makeRFs(washedDat[['expDat']][,rownames(stTrain)], stTrain, geneLists, dLevel="bottomL",
                     nTrees=nTrees)
    
    # classify held out data
    stVal<-ttList[['stVal']]
    ansVal=sc_classify(testRFs, washedDat[['expDat']][,rownames(stVal)], geneLists)
    
    
    cat("Assessing bottomLevel\n")
    assessed<-sc_classAssess(ansVal, stVal, classLevels="bottomL", dLevelSID=dLevelSID, resolution=0.01)
    classifiers[['bottomLevel']] <- list(classifiers=testRFs, classRes=ansVal, PRs=assessed, stVal = stVal)
  }
  
  
  return(classifiers) 
}  


sc_classAssess<-function# make ROCs for each classifier
(ansVal,# matrix of classification scores, rows = classifiers, columns = samples, colnames=sampleids
 stVal, # sampletable
 classLevels=classLevels,
 dLevelSID=dLevelSID,
 resolution=0.005 # increment at which to evalutate classification
){
  allROCs<-list();
  evalAll<-matrix(0, nrow=nrow(ansVal),ncol=2)
  classifications<-rownames(ansVal)
  for(i in 1:length(classifications)){
    classification<-classifications[i];
    tmpROC <- cn_eval(ansVal[classification,],
                      stVal,
                      classLevels,
                      classification,threshs=seq(0,1, by=resolution), dLevelSID=dLevelSID);
    allROCs[[classification]]<-tmpROC;
  }
  allROCs;
}

cn_eval<-function# return a data frame of the number of TP, FN, FP, and TN, and pval cutoff
(vect, # named vector
 sampTab,
 dLevel, # description level)
 classification,
 threshs=seq(0,1,by=0.05), # pval cutoffs
 dLevelSID=dLevelSID
){
  ans<-matrix(0,nrow=length(threshs), ncol=7);
  for(i in seq(length(threshs))){
    thresh<-threshs[i];
    ans[i,1:4]<-cn_clPerf(vect, sampTab, dLevel, classification, thresh, dLevelSID=dLevelSID);
  }
  ans[,5]<-threshs;
  colnames(ans)<-c("TP", "FN", "FP", "TN", "thresh","FPR", "TPR");
  TPR<-ans[,'TP']/(ans[,'TP']+ans[,'FN']);
  FPR<-ans[,'FP']/(ans[,'TN']+ans[,'FP']);
  ans[,'TPR']<-TPR;
  ans[,'FPR']<-FPR;
  ans;
}
cn_clPerf<-function # assumes rownames(sampTab) == sampTab identifier used as colname for vect
(vect,
 sampTab,
 dLevel,
 classification, # actual classification
 thresh,
 dLevelSID=dLevelSID){
  TP<-0;
  FN<-0;
  FP<-0;
  TN<-0;
  sampIDs<-names(vect);  
  classes<-as.vector(sampTab[sampIDs,dLevel]);
  
  ###actualPos<-as.vector(sampTab[sampTab[,dLevel]==classification,]$sample_id);#which(classes==classification));
  actualPos<-as.vector(sampTab[sampTab[,dLevel]==classification,dLevelSID])
  actualNeg<-setdiff(sampIDs, actualPos);
  
  calledPos<-names(which(vect>thresh));
  calledNeg<-names(which(vect<=thresh));
  
  TP <- length(intersect(actualPos, calledPos));
  FP <- length(intersect(actualNeg, calledPos));
  FN <- length(actualPos)-TP;
  TN <- length(actualNeg)-FP;
  c(TP, FN, FP, TN);  
}


gnrAll<-function(
  expDat,
  cellLabels){

  myPatternG<-sc_sampR_to_pattern(as.character(cellLabels))
  specificSets<-lapply(myPatternG, sc_testPattern, expDat=expDat)
  cat("Done testing\n")
  specificSets
}

#' getClassGenes
#'
#' extract genes for training classifier
#' @param diffRes a df with cval, holm, rownames=genes
#' @param topX number of genes to select
#' @param bottom boolean if ture use the top x genes with - cvals
#'
#' @return vector of genes
#'
#' @export
getClassGenes<-function(
 diffRes,
 topX=25,
 bottom=TRUE)
 {
   #exclude NAs
   xi<-which(!is.na(diffRes$cval))
   diffRes<-diffRes[xi,]   
   diffRes<-diffRes[order(diffRes$cval, decreasing=TRUE),]
   ans<-rownames(diffRes[1:topX,])
   if(bottom){
     ans<-append(ans, rownames( diffRes[nrow(diffRes) - ((topX-1):0),]))
   }
   ans
 }