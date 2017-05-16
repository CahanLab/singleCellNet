# scCellNet
# (C) Patrick Cahan 2012-2017

#' make classifiers
#'
#' Make one Random Forest classification per unique dLevel
#' @param expTrain expression matrix
#' @param stTrain sample table
#' @param list of dLevel=>vector of genes
#' @param dLevel stTrain column name to split on
#'
#' @return list of classifiers
makeRFs<-function#
(expTrain,
 stTrain,
 geneLists, 
 dLevel='description1'){
  ans<-list();
  cnames<-names(geneLists);
  cnames<-intersect(cnames, unique(as.vector(stTrain[,dLevel])));
  for(cname in cnames){
    cat(cname,"\n");
    xgenes<-intersect(rownames(expTrain), geneLists[[cname]]);  
    cat(cname ,":", length(xgenes),"\n");
    ans[[cname]]<-makeClassifier(expTrain[xgenes,],cname,stTrain, dLevel=dLevel);
  }
  ans;
}

#' make a single RF classifier
#'
#' uses R's randomForest package
#' @param trainScores numeric matrix of predictors
#' @param ctt what cell type/tissues is this?
#' @param stTrain sample table
#'
#' @return randomForest classifier
makeClassifier<-function 
(trainScores, 
 ctt,
 stTrain,
 dLevel='description1'){
  resp<-as.vector(stTrain[,dLevel]);
  xi<-which(resp!=ctt);
  resp[xi]<-'other';  
  myClass<-randomForest(t(trainScores), as.factor(resp), ntree=2e3);
  myClass;
}

#' classify data
#'
#' run binary classifiers on input expression matrix
#' @param classList result of running CN3_makeClassifier
#' @param expDat  expression data matrix
#' @param cttComms list of genes that were used to train each classifer cttComms[[className]]<-c(gene1, gene2, ...)
#'
#' @return classification matrix nrow=length(classList) ncol=ncol(expDat)
#'
sc_classify<-function
(classList, 
 expDat, 
 cttComms 
){
  ctts<-names(classList);
  ans<-matrix(0, nrow=length(ctts), ncol=ncol(expDat));
  rownames(ans)<-ctts;
  for(ctt in ctts){
    #cat(ctt,"\n")
    ##    xgenes<-cttComms[[ctt]][[1]];
    xgenes<-cttComms[[ctt]];
    ### 06-05-16
    xgenes<-intersect(xgenes, rownames(expDat));
    myProbs<-predict(classList[[ctt]], t(expDat[xgenes,]), type='prob');
    ans[ctt,]<-myProbs[,ctt];
  }
  colnames(ans)<-colnames(expDat);
  ans;
  # classification matrix
}

#' split data into train vs test
#'
#' Splits a sample table into 2 sample tables roughly by prop in which no samples with sampe exp_id are in both the test and train
#' @param sampTab sample table must have exp_id column
#' @param prop proportion of samples to use as training
#' @param dLevel column name for 
#' @param dLevelStudy column name to indicate experiment or study id
#'
#' @return list of stTrain stVal
samp_for_class<-function# return a sampTab for training a test classifier
(sampTab,
  prop=0.5,
  dLevel="description1",
  dLevelStudy="exp_id"){

  ctts<-unique(as.vector(sampTab[,dLevel]));
  stTrain<-data.frame();
  stVal<-data.frame();
  for(ctt in ctts){
    stTmp<-sampTab[sampTab[,dLevel]==ctt,];
    stTrain<-rbind(stTrain, subSamp_for_class(stTmp,prop, dLevelStudy));

    idsval<-setdiff( rownames(stTmp), rownames(stTrain) );
    stVal<-rbind(stVal, stTmp[idsval,]);
  }
  list(stTrain=stTrain, stVal=stVal);
}

#' select from sample table
#'
#' that's it. helper function
#' @param sampTab sample table
#' @param prop fraction of samples to select
#' @param dLevelStudy column name to indicate experiment or study id
#'
#' @return stTrain
subSamp_for_class<-function
(sampTab,
  prop=0.5,
  dLevelStudy="exp_id"){

  expIDcounts<-sort(table(sampTab[,dLevelStudy]));
  expIDs<-names(expIDcounts);
  total<-sum(expIDcounts);

  runningTotal<-0;
  i<-1;
  xi<-floor( length(expIDcounts)/2 );
  while(i<=length(expIDcounts)){
    runningTotal<-sum(expIDcounts[1:i]);
    if(runningTotal/total > prop){
      xi<- i-1;
      break;
    }
    i<-i+1;
  }
  expIDs<-expIDs[1:xi];
  
  stTrain<-data.frame();
  for(expID in expIDs){
    stTrain<-rbind(stTrain, sampTab[sampTab[,dLevelStudy]==expID,]);
  }
  stTrain;
}




