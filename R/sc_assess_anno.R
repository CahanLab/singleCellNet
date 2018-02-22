#' assess the classifiers performance based on validation data
#' @param stDat   sample table
#' @param expDat  normalized expression matrix
#' @param nimCells  the minimal number of cells one would like to have in each cell type
#' @param dThresh  detection threshold
#' @param propTrain  the proportion of the training data desire
#' @param nRand  the number of random sample one wants to generate 
#' @param nTrees number of branches one would like to build on the random forest classifier 


sc_classAssess<-function
(stDat,
 washedDat,
 dLevel = "description1",
 dLevelSID="sample_name",
 minCells = 40,
 dThresh = 0, 
 propTrain=0.25,
 nRand = 50,
 nTrees=2000,
 resolution=0.005){
  
  #extract sample table and gene expression matrix of groups that are above minCells  
  goodGrps<-names(which(table(stDat[,dLevel])>=minCells))
  stTmp<-data.frame()
  for(ggood in goodGrps){
    stTmp<-rbind(stTmp, stDat[stDat[,dLevel]==ggood,])
  }
  dim(stTmp)
  expDat_good <- washedDat$expDat[,rownames(stTmp)]

  # split dataset into training and test data
  ttList<-divide_sampTab(stTmp, propTrain, dLevel = dLevel)
  stTrain<-ttList[['stTrain']]
  expTrain <- expDat_good[,row.names(stTrain)]

    #building classifiers on varGenes
  #can build in options to select different methods
  cat("Calculating varGenes...\n")
  varGenes <- findVarGenes(expDat_good, washedDat$geneStats)

  cellgrps<-stTrain[,dLevel]
  names(cellgrps)<-rownames(stTrain)

  cat("Making classifiers...\n")
  testRFs<-sc_makeClassifier(expTrain, genes=varGenes, groups=cellgrps, nRand=nRand, ntrees=nTrees)

  cat("classify held out data...\n")
  stVal<-ttList[['stVal']]
  ct_scores=rf_classPredict(testRFs, expDat_good[,row.names(stVal)])

  #cat("obtaining some ROCs...\n")
  #assessed<-sc_ROCs(ct_scores, stVal, dLevel=dLevel, dLevelSID=dLevelSID, threshs=seq(0,1,by=resolution, nRand = nRand))

}  

#' Assess classifiers based on validation data
#'
#' Assess classifiers based on validation data
#' @param ct_scores matrix of classification scores, rows = classifiers, columns = samples, colnames=sampleids
#' @param stVal sample table
#' @param classLevels column name of stVal to use as ground truth to assess classifiers
#' @param threshs seq of pval cutoffs
#' @param dLevelSID column to indicate sample id
#'
#' @return list of data frames with threshold, sens, precision
#' @export 

sc_ROCs<-function# make ROCs for each classifier, pairwise comparison
(ct_scores,# matrix of classification scores, rows = classifiers, columns = samples, colnames=sampleids
 stVal, # sampletable
 dLevel="description1",
 dLevelSID="sample_name",
 threshs=seq(0,1,by=0.005),
 nRand = 50 # increment at which to evalutate classification
){

  sampIDs<-colnames(ct_scores)

  #make a new stVal with the rand cells added
  tmp <- as.data.frame(matrix("Rand", nrow = nRand, ncol=(ncol(stVal))))
  colnames(tmp) <- colnames(stVal)
  tmp$sample_name <- sampIDs[(length(sampIDs) - nRand + 1):length(sampIDs)]
  rownames(tmp) <- tmp[,dLevelSID]
  stVal_tmp <- rbind(stVal, tmp)
  
  ans<-matrix(0,nrow=length(threshs), ncol=7);



  for(i in seq(length(threshs))){
    thresh<-threshs[i];
    ans[i,1:4]<-sc_Perf(ct_scores, stVal_tmp, dLevel, thresh, dLevelSID=dLevelSID);
  }
  ans[,5]<-threshs;
  colnames(ans)<-c("TP", "FN", "FP", "TN", "thresh","FPR", "TPR");
  TPR<-ans[,'TP']/(ans[,'TP']+ans[,'FN']);
  FPR<-ans[,'FP']/(ans[,'TN']+ans[,'FP']);
  ans[,'TPR']<-TPR;
  ans[,'FPR']<-FPR;
  ans;

}


#' determine performance of classification at given threshold
#'
#' determine performance of classification at given threshold
#' @param ct_scores matrix of classification scores, rows = classifiers, columns = samples, colnames=sampleids
#' @param sampTab sample table
#' @param dLevel colname
#' @param classification actual classification
#' @param thresh threshold above which to make a call
#' @param dLevelSID column to indicate sample id
#'
#' @return vector of TP FN FP TN 
sc_Perf<-function # assumes rownames(sampTab) == sampTab identifier used as colname for vect
(ct_scores,
 sampTab,
 dLevel,
 classification, # actual classification
 thresh,
 dLevelSID="sample_name"){
  TP<-0;
  FN<-0;
  FP<-0;
  TN<-0;

  sampIDs<-colnames(ct_scores);  
  
  ###actualPos<-as.vector(sampTab[sampTab[,dLevel]==classification,]$sample_id);#which(classes==classification));
  actualPos<-as.vector(sampTab[sampTab[,dLevel],dLevelSID])
  actualNeg<-setdiff(sampIDs, actualPos);
  calledPos<-vector()
  calledNeg<-vector() 
  
  classifications<-rownames(ct_scores);

  for (classification in classifications){
    vect <- ct_scores[rownames(ct_scores) == classification,]
    calledPos <- c(calledPos, names(which(vect>thresh)))
    calledNeg<-c(calledNeg, names(which(vect<=thresh)))

  }
  
  TP <- length(intersect(actualPos, calledPos));
  FP <- length(intersect(actualNeg, calledPos));
  FN <- length(actualPos)-TP;
  TN <- length(actualNeg)-FP;
  c(TP, FN, FP, TN);  
}


#' Plot results of sc_classAssess
#'
#' Plot one precision recall curve per CT
#' @param assessed result of runnung cn_classAssess
#'
#' @return ggplot pbject
#'
#' @examples
#' testAssTues<-cn_splitMakeAssess(stTrain, expTrain, ctGRNs, prop=.5)
#' plot_class_PRs(testAssTues$ROCs)
#'
#' @export
plot_class_PRs<-function
(assessed
  ){
  ctts<-names(assessed);
  df<-data.frame();
  for(ctt in ctts){
    tmp<-assessed[[ctt]];
    tmp<-cbind(tmp, ctype=ctt);
    df<-rbind(df, tmp);
  }

  prsAll<-transform(df, TP = as.numeric(as.character(TP)), 
    TN = as.numeric(as.character(TN)), 
    FN = as.numeric(as.character(FN)), 
    FP = as.numeric(as.character(FP)));

    precfunc<-function(df){
      ans<-vector();
      for(i in 1:nrow(df)){
        ans<-append(ans, df[i,"TP"]/(df[i,"TP"]+df[i,"FP"]));
      }
      ans;
    }

    sensfunc<-function(df){
      ans<-vector();
      for(i in 1:nrow(df)){
        ans<-append(ans, df[i,"TP"]/(df[i,"TP"]+df[i,"FN"]));
      }
      ans;
    }

  precs<-precfunc(prsAll)
  sens<-sensfunc(prsAll)
  prsAll2<-cbind(prsAll, data.frame(recall=sens, precision=precs));

  ggplot(data=prsAll2, aes(x=as.numeric(as.vector(recall)), y=as.numeric(as.vector(precision)))) + geom_point(size = .5, alpha=.5) +  geom_path(size=.5, alpha=.75) +
  theme_bw() + xlab("Recall") + ylab("Precision") + facet_wrap( ~ ctype, ncol=4) +
  theme(axis.text = element_text(size=5)) + ggtitle("Classification performance")
}

findVarGenes<-function(expNorm, geneStats,zThresh=2){
  allGenes<-rownames(expNorm)
  sg<-binGenesAlpha(geneStats)
  zscs<-rep(0, nrow(sg));
  names(zscs)<-rownames(sg);
  bbins<-unique(sg$bin);
  for(bbin in bbins){
    xx<-sg[sg$bin==bbin,];
    tmpZ<-scale(xx$fano);
    zscs[ rownames(xx) ]<-tmpZ[,1];
  }

  zByAlpha<-names(which(zscs>zThresh))
  
  # by overall mean
  sg<-binGenes(geneStats)
  zscsM<-rep(0, nrow(sg));
  names(zscsM)<-rownames(sg);
  bbins<-unique(sg$bin);
  for(bbin in bbins){
    xx<-sg[sg$bin==bbin,];
    tmpZ<-scale(xx$fano);
    zscsM[ rownames(xx) ]<-tmpZ[,1];
  }
  zByMean<-names(which(zscsM>zThresh))
  union(zByMean, zByAlpha)
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

 binGenesAlpha<-function
(geneStats,
 nbins=20){
  max<-max(geneStats$alpha);
  min<-min(geneStats$alpha);

  binGroup<-rep(nbins, length=nrow(geneStats));
  names(binGroup)<-rownames(geneStats);
  rrange<-max-min;
  inc<-rrange/nbins
  borders<-seq(inc, max, by=inc)
  for(i in length(borders):1){
    xnames<-rownames(geneStats[which(geneStats$alpha<=borders[i]),]);
    binGroup[xnames]<-i;
  }
  cbind(geneStats, bin=binGroup);
}

### bins genes into x groups based on overallmean
binGenes<-function
(geneStats,
 nbins=20,
 meanType="overall_mean"){

##  max<-max(geneStats$overall_mean);
##  min<-min(geneStats$overall_mean);

  max<-max(geneStats[,meanType])
  min<-min(geneStats[,meanType])

  binGroup<-rep(nbins, length=nrow(geneStats));
  names(binGroup)<-rownames(geneStats);
  rrange<-max-min;
  inc<-rrange/nbins
  borders<-seq(inc, max, by=inc)
  for(i in length(borders):1){
    #xnames<-rownames(geneStats[which(geneStats$overall_mean<=borders[i]),]);
    xnames<-rownames(geneStats[which(geneStats[,meanType]<=borders[i]),]);
    binGroup[xnames]<-i;
  }
  cbind(geneStats, bin=binGroup);
}

#' compute AUPCR
#'
#' compute AUPCR
#' @param perfDF a df
#' @param precisionCol "Precision"
#' @param recallCol "Recall"
#' @param predCol "Predictions"
#'
#' @return auprs
cn_computeAUCPR<-function
(perfDF,
 precisionCol="Precision",
 recallCol="Recall",
 predCol="Predictions"
){
  
  ### Notes: starts at top left, and accumulates to max
  str<-(nrow(perfDF)-1);
  stp<-2;
  
  # sometimes you can still get NA in 
  areas<-rep(0, nrow(perfDF));
  
  pts<-seq(str,stp,-1);
  for(i in pts){
    a<-(i+1);
    cc<-(i-1);
    ptA<-c(perfDF[a,recallCol], perfDF[a,precisionCol]);
    ptB<-c(perfDF[i,recallCol], perfDF[i,precisionCol]);
    ptC<-c(perfDF[cc,recallCol], perfDF[cc,precisionCol]);
    tmpArea<-cn_rectArea(ptA, ptB, ptC);
    if(is.nan(tmpArea)){
      #cat(perfDF[i,]$Score,"\n");
      tmpArea<-0;
    }
    areas[i]<-areas[(i+1)]+tmpArea;
  }
  
  # far right point
  a<-2;
  b<-1;
  cc<-1;
  ptA<-c(perfDF[a,recallCol], perfDF[a,precisionCol]);
  ptB<-c(perfDF[i,recallCol], perfDF[i,precisionCol]);
  ptC<-c(perfDF[cc,recallCol], perfDF[cc,precisionCol]);
  areas[b]<-areas[2]+cn_rectArea(ptA, ptB, ptC);
  
  # far left point
  a<-nrow(perfDF);
  b<-nrow(perfDF);
  cc<-(b-1);
  ptA<-c(perfDF[a,recallCol], perfDF[a,precisionCol]);
  ptB<-c(perfDF[i,recallCol], perfDF[i,precisionCol]);
  ptC<-c(perfDF[cc,recallCol], perfDF[cc,precisionCol]);
  areas[b]<-cn_rectArea(ptA, ptB, ptC);
  areas;
}

#' compute area of rect given by 
#'
#' compute area of rect given by 
#' @param ptA point a
#' @param ptB point b
#' @param ptC point C
#'
#' @return area
cn_rectArea<-function
(ptA,
 ptB,
 ptC
){
  xRight <- ptB[1]+( (ptC[1]-ptB[1])/2);
  xLeft  <- ptA[1]+( (ptB[1]-ptA[1])/2);
  width<-abs(xRight - xLeft);
  width*ptB[2];
}
