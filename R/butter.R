# singleCellNet
# (C) Patrick Cahan 2012-2017


extractImport<-function
(classifiers # from makeRFs
){
  classNames<-names(classifiers)
  preds<-rownames(classifiers[[1]]$importance)
  ans<-matrix(0, nrow=length(preds), ncol=length(classNames))
  rownames(ans)<-preds
  colnames(ans)<-classNames
  for(cname in classNames){
    ans[,cname]<-classifiers[[cname]]$importance
  }
  ans
}


#' @export
butter_classGenes<-function
(washed, # from wash()
 sampTab, # has group col, likely from steam
 classifiers, # from pipe_butter
 ngenes=10 # number of genes per class
 ){

  diffGenesMat<-enrDiff(washed[['expDat']], sampTab) # returns a matrix, one col per cluster
  classList<-classifiers$classifiers
  importGenes<-extractImport(classList)
  diffGenesMat<-diffGenesMat[rownames(importGenes),] * (1+importGenes)
  classNames<-names(classList)
  ans<-list()
  for(cname in classNames){
#    ans[[cname]]<-rownames(diffGenesMat)[order(diffGenesMat[,cname], decreasing=TRUE)][1:ngenes]

    tmpans<-diffGenesMat[order(diffGenesMat[,cname], decreasing=TRUE),cname]
    xi<-which(tmpans<=0)[1]
    tmpThresh<-min(ngenes, (xi-1))
    ans[[cname]]<-names(tmpans[1:tmpThresh])
  }
  ans

}


#' classify
#'
#' classify, wrapper to sc_classify
#' @param washedDat result of wash()
#' @param classList result of pipeButter()
#'
#' @return matrix of classifier results
#'
#' @export
butter_classify<-function
(washedDat,
 classList # result of pipeButter
 ){

    # fill out missing genes
    predictors<-classList[['predictors']][[1]]
    genes<-rownames(washedDat[['expDat']])
    mGenes<-setdiff(predictors, genes)
    cat(mGenes,"\n")
    matToAdd<-matrix(0, nrow=length(mGenes), ncol=ncol(washedDat[['expDat']]))
    rownames(matToAdd)<-mGenes
    colnames(matToAdd)<-colnames(washedDat[['expDat']])
    washedDat[['expDat']]<-rbind(washedDat[['expDat']], matToAdd)

    # !!! add check that wash is done the same way
    sc_classify(classList[['classifiers']], washedDat[['expDat']], classList[['predictors']])
}

#' make classifiers
#'
#' Make one Random Forest classification per unique dLevel
#' @param expTrain expression matrix
#' @param stTrain sample table
#' @param list of dLevel=>vector of genes
#' @param dLevel stTrain column name to split on
#'
#' @return list of classifiers
#'
#' @export
makeRFs<-function#
(expTrain,
 stTrain,
 geneLists, 
 dLevel='description1',
 nTrees=2e2){
  ans<-list();
  cnames<-names(geneLists);
  cnames<-intersect(cnames, unique(as.vector(stTrain[,dLevel])));
  for(cname in cnames){
   ### cat(cname,"\n");
    xgenes<-intersect(rownames(expTrain), geneLists[[cname]]);  
    ### cat(cname ,":", length(xgenes),"\n");
    ans[[cname]]<-makeClassifier(expTrain[xgenes,],cname,stTrain, dLevel=dLevel,nTrees)
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
#'
#' @export
makeClassifier<-function 
(trainScores, 
 ctt,
 stTrain,
 nTrees=2e2,
 dLevel='description1'){
  resp<-as.vector(stTrain[,dLevel]);
  xi<-which(resp!=ctt);
  resp[xi]<-'other';  
  myClass<-randomForest(t(trainScores), as.factor(resp), ntree=nTrees)
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
#' @export
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

#' divide sample table
#'
#' divide sample table
#' @param sampTab sampTab
#' @param prop prop
#' @param dLevel dLevel
#'
#' @return list of stTrain stVal
#'
#' @export
divide_sampTab<-function
(sampTab,
  prop=0.5,
  dLevel="prefix"
){

  ctts<-unique(as.vector(sampTab[,dLevel]));
  stTrain<-data.frame();
  stVal<-data.frame();

  for(ctt in ctts){
    stTmp<-sampTab[sampTab[,dLevel]==ctt,];
    xCount<-ceiling(prop * nrow(stTmp))
    xnames<-sample(rownames(stTmp), xCount)

    stTrain<-rbind(stTrain, stTmp[xnames,])

    idsval<-setdiff( rownames(stTmp), rownames(stTrain) );
    stVal<-rbind(stVal, stTmp[idsval,]);
  }
  list(stTrain=stTrain, stVal=stVal);
}


#' @export
easy_assess<-function
(classRes,
 sampTab)
{
   afunct<-function(vector, len){
     xi<-which.max(vector)
     others<-setdiff(1:len, xi)
     vector[xi]-sum(vector[others])
   }
   ans<-apply(classRes, 2, afunct, nrow(classRes))
   cbind(sampTab, classDiff=ans)
 }

# splits data
# makes classifiers
# apply to held out data
#' @export
prebutter<-function
(predictors,
  expDat,
  stDat, # with category partOn
 propTrain=0.25,
 partOn="group",
 nTrees=200)
{
  geneLists<-list()
  cts<-as.vector(unique(stDat[,partOn]))
  for(ct in cts){
    geneLists[[ct]]<-predictors
  }

  # split into training and test data
  ttList<-divide_sampTab(stDat, propTrain, partOn)
  stTrain<-ttList[['stTrain']]

  # make RFs
  myRFs<-makeRFs(expDat[predictors,rownames(stTrain)], stTrain, geneLists, dLevel=partOn,
    nTrees=nTrees)
  
  # classify held out data
  stVal<-ttList[['stVal']]

  list(classRes=sc_classify(myRFs, expDat[predictors,rownames(stVal)], geneLists), stVal=stVal)
}




















