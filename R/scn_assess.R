#' assess the classifiers performance based on validation data
#' @param stDat   sample table
#' @param expDat  normalized expression matrix
#' @param nimCells  the minimal number of cells one would like to have in each cell type
#' @param dThresh  detection threshold
#' @param propTrain  the proportion of the training data desire
#' @param nRand  the number of random sample one wants to generate 
#' @param nTrees number of branches one would like to build on the random forest classifier 
#' @export 

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

  cellgrps<-as.character(stTrain[,dLevel])
  names(cellgrps)<-rownames(stTrain)

  cat("Making classifiers...\n")
  testRFs<-sc_makeClassifier(expTrain, genes=varGenes, groups=cellgrps, nRand=nRand, ntrees=nTrees)

  cat("classify held out data...\n")
  stVal<-ttList[['stVal']]
  ct_scores=rf_classPredict(testRFs, expDat_good[,row.names(stVal)])

  assessed <- list(ct_scores = ct_scores, stVal = stVal, stTrain = stTrain)
  
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
#' @return list of data frames with threshold
#' @export 

sc_Accu<-function #calculate the accuracy of the data at each given classification threshold
(ct_scores,# matrix of classification scores, rows = classifiers, columns = samples, colnames=sampleids
 stVal, # sampletable
 dLevel="description1",
 dLevelSID="sample_name",
 threshs=seq(0,1,by=0.005),
 nRand = 50 # increment at which to evalutate classification
){

  sampIDs<-colnames(ct_scores)

  #make a new stVal with the rand cells added
  colnames(tmp) <- colnames(stVal)
  tmp[,dLevelSID] <- sampIDs[(length(sampIDs) - nRand + 1):length(sampIDs)]
  rownames(tmp) <- tmp[,dLevelSID]
  stVal_tmp <- rbind(stVal, tmp)
  
  ans<-matrix(0,nrow=length(threshs), ncol=nrow(stVal_tmp));
  rownames(ans) <- threshs; #rownames is the classification threshold
  colnames(ans) <- sampIDs

  for (i in 1: length(sampIDs)){
    sampID <- sampIDs[i]
    vect <- ct_scores[,colnames(ct_scores) == sampID]
    classification <- as.vector(stVal_tmp[sampID, dLevel])

    for (j in 1: length(threshs)){  
      thresh <- threshs[j]
      ans[j,sampID] <- sc_classThreshold(vect, classification, thresh)
    }

  }
  
  ans;

}


#' determine performance of classification at given class threshold
#'
#' determine performance of classification at given threshold
#' @param vect, classifcation score for each given cell
#' @param classification actual classification
#' @param thresh threshold above which to make a call
#'
#' @return accuracy
#' @export 

sc_classThreshold<-function # assumes rownames(sampTab) == sampTab identifier used as colname for vect
(vect,
 classification, #actual classification
 thresh)
{

  TP<-0;
  FN<-0;
  FP<-0;
  TN<-0;
  
  calledPos <- names(which(vect>thresh))
  calledNeg <- names(which(vect<=thresh))

  if (classification %in%  calledPos){
      TP = 1
      FN = 0
      FP = length(calledPos) - 1
      TN = length(calledNeg)
    } else {
      TP = 0
      FN = 1
      FP = length(calledPos)
      TN = length(calledNeg) -1
    }

  Accu <- (TP + TN)/(TP + TN + FP + FN)

}

plot_accu <- function
(accu_dat, #data from sc_Accu output
 threshs){

newDat <- cbind(threshs, accu_dat)
newDat <- as.data.frame(newDat)

#melt the data by the first colomn
dat.m = melt(newDat, id = "threshs")

col=colorRampPalette(rev(brewer.pal(n = 11,name = "RdYlBu")))(length(threshs))

ggplot(data = dat.m,aes(x=as.character(threshs), y=value)) + geom_boxplot(fill = col) +labs(title="Plot of Accuracy per Classification Threshold",x="Classification Thresholds", y = "Accuracy")

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


#' compute AUPCR
#'
#' compute AUPCR
#' @param perfDF a df
#' @param precisionCol "Precision"
#' @param recallCol "Recall"
#' @param predCol "Predictions"
#'
#' @return auprs
#' @export 
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
#' @export 
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

#assessment
#' Assess classifiers based on validation data
#'
#' Assess classifiers based on validation data
#' @param ct_scores matrix of classification scores, rows = classifiers, columns = samples, colnames=sampleids
#' @param stVal sample table
#' @param classLevels column name of stVal to use as ground truth to assess classifiers
#' @param resolution increment at which to evalutate classification
#' @param dLevelSID column to indicate sample id
#'
#' @return list of data frames with threshold, sens, precision
#' @export 
cn_classAssess<-function# make ROCs for each classifier
(ct_scores,# matrix of classification scores, rows = classifiers, columns = samples, colnames=sampleids
 stVal, # sampletable
 classLevels="description2",
 dLevelSID="sample_id",
 resolution=0.005 # increment at which to evalutate classification
){
  allROCs<-list();
  evalAll<-matrix(0, nrow=nrow(ct_scores),ncol=2);
  classifications<-rownames(ct_scores);
  rownames(stVal)<-as.vector(stVal[,dLevelSID]);
  i<-1;
  for(xname in classifications){
    classification<-classifications[i];
    tmpROC <- cn_eval(ct_scores[xname,],
                           stVal,
                           classLevels,
                           xname,threshs=seq(0,1, by=resolution), dLevelSID=dLevelSID);
    allROCs[[xname]]<-tmpROC;
    i<-1+i;
  }
  allROCs;
}

#' return the sens at given FPR
#'
#' return the sens at given FPR
#' @param rocRes result of running cn_eval
#' @param what is sens at this FPR?
#'
#' @return what is sens at selected FPR?
cn_findSensAt<-function
(rocRes,
 at=0.05){
  xi<-which(rocRes[,"FPR"]<=at);
  tmp<-rocRes[xi,];
  cat(nrow(tmp),"\n")
  as.vector(tmp[which.max(tmp[,"FPR"]),'TPR']);
}

#' run cn_clPerf across thresholds
#'
#' run cn_clPerf across thresholds
#' @param vect named vector
#' @param dLevel description level
#' @param classification classification matrix
#' @param threshs seq of pval cutoffs
#' @param dLevelSID column to indicate sample id
#'
#' @return return a data frame of the number of TP, FN, FP, and TN, and pval cutoff
cn_eval<-function# return a data frame of the number of TP, FN, FP, and TN, and pval cutoff
(vect, # named vector
 sampTab,
 dLevel, # description level)
 classification,
 threshs=seq(0,1,by=0.05), # pval cutoffs
 dLevelSID="sample_id"
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

#' determine performance of classification at given threshold
#'
#' determine performance of classification at given threshold
#' @param vect vector of values
#' @param sampTab sample table
#' @param dLevel colname
#' @param classification actual classification
#' @param thresh threshold above which to make a call
#' @param dLevelSID column to indicate sample id
#'
#' @return vector of TP FN FP TN 
cn_clPerf<-function # assumes rownames(sampTab) == sampTab identifier used as colname for vect
(vect,
 sampTab,
 dLevel,
 classification, # actual classification
 thresh,
 dLevelSID="sample_id"){
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

#' @param ct_score cnRes score
#' @param stQuery query sample table
#' @param nRand number of rand used in rand
#' @param dLevelSID columns used to split the samples
#' @export
makeSampleTable <- function(ct_scores, stQuery, nRand, dLevelSID){
  
  sampIDs<-colnames(ct_scores)
  tmp <- as.data.frame(matrix("rand", nrow = nRand, ncol=(ncol(stQuery))))
  colnames(tmp) <- colnames(stQuery)
  tmp[,dLevelSID] <- sampIDs[(length(sampIDs) - nRand + 1):length(sampIDs)]
  rownames(tmp) <- tmp[, dLevelSID]
  stVal_tmp <- rbind(stQuery, tmp)
  rownames(stVal_tmp) <- stVal_tmp[,dLevelSID]
  return(stVal_tmp)
}

#' @export
assess_comm <- function(ct_scores, #matrix of classification scores, rows = classifiers, columns = samples, colnames=sampleids || where query cells is in the training
         stTrain, #sample table where cells in query are in the training 
         stQuery,
         resolution = 0.005,# increment at which to evalutate classification
         nRand = 50,
         dLevelSID = "sample_name",
         classTrain = "cell_ontology_class",
         classQuery = "description2", #query data
         AUCmethod = "trapezoid"){
  
  shared_cell_type <- intersect(unique(stTrain[,classTrain]), unique(stQuery[,classQuery]))
  stVal_com <- stQuery[which(stQuery[,classQuery] %in% shared_cell_type),]
  
  if(nRand > 0){
    tmp <- as.data.frame(matrix("rand", nrow = nRand, ncol=(ncol(stVal_com))))
    colnames(tmp) <- colnames(stVal_com)
    tmp[,dLevelSID] <- colnames(ct_scores)[(ncol(ct_scores) - nRand + 1):ncol(ct_scores)]
    rownames(tmp) <- tmp[,dLevelSID]
    stVal_com <- rbind(stVal_com, tmp)
  }
  
  cells_sub <- as.character(stVal_com[,dLevelSID])
  
  #subsetting the ct_scores where the cells' true identity is within the range of the classifiers
  ct_score_com <- ct_scores[, cells_sub]
  report <- list()
  ct_scores_t <- t(ct_score_com)
  true_label <- as.character(stVal_com[, classQuery])
  #multiLogLoss
  names(true_label) <- rownames(ct_scores_t)
  if (is.matrix(true_label) == FALSE) {
    y_true <- model.matrix(~ 0 + ., data.frame(as.character(true_label)))
  }
  eps <- 1e-15
  y_pred <- pmax(pmin(ct_scores_t, 1 - eps), eps)
  multiLogLoss <- (-1 / nrow(ct_scores_t)) * sum(t(y_true)%*% log(y_pred)) #want columns to be the cell types for y_pred
  #cohen's kappa, accuracy
  pred_label <- c()
  pred_label <- colnames(ct_scores_t)[max.col(ct_scores_t,ties.method="random")]
  
  cm = as.matrix(table(Actual = true_label, Predicted = pred_label))
  
  #in case of misclassfication where there are classifiers that are not used
  if(length(setdiff(unique(true_label), unique(pred_label))) != 0){
    misCol <- setdiff(unique(true_label), unique(pred_label))
    for(i in 1:length(misCol)){
      cm <- cbind(cm, rep(0, nrow(cm)))
    }
    colnames(cm)[(ncol(cm) - length(misCol) +1) : ncol(cm)] <- misCol
  }
  
  if(length(setdiff(unique(pred_label), unique(true_label))) != 0){
    misRow <- setdiff(unique(pred_label), unique(true_label))
    for(i in 1:length(misRow)){
      cm <- rbind(cm, rep(0, ncol(cm)))
    }
    rownames(cm)[(nrow(cm) - length(misRow) +1) : nrow(cm)] <- misRow
  }
  
  cm <- cm[,colnames(cm)[match(rownames(cm),colnames(cm))]]  
  
  
  #sort table names accordigly
  
  n = sum(cm) # number of instances
  nc = nrow(cm) # number of classes
  diag = diag(cm) # number of correctly classified instances per class 
  rowsums = apply(cm, 1, sum) # number of instances per class
  colsums = apply(cm, 2, sum) # number of predictions per class 
  p = rowsums / n # distribution of instances over the actual classes
  q = colsums / n # distribution of instances over the predicted classes
  expAccuracy = sum(p*q)
  accuracy = sum(diag) / n 
  
  #PR
  confusionMatrix <- cn_classAssess(ct_score_com, stVal_com, classLevels= classQuery, dLevelSID=dLevelSID, resolution=resolution)
  PR_ROC <- cal_class_PRs(confusionMatrix)
  nonNA_PR <-PR_ROC[which(!is.nan(PR_ROC$recall)),]
  nonNA_PR[which((nonNA_PR$TP == 0 & nonNA_PR$FP ==0)), "precision"] <- 1
  
  w <- c()
  areas <- c()
  for(i in 1: length(unique(nonNA_PR$ctype))){
    tmp <- nonNA_PR[which(nonNA_PR$ctype %in% unique(nonNA_PR$ctype)[i]),]
    area <- DescTools::AUC(tmp$recall, tmp$precision, method = AUCmethod)
    areas <- c(areas,area[1])
    w <- c(w, sum(stVal_com[,classQuery] %in% unique(nonNA_PR$ctype)[i])/nrow(stVal_com))
  }
  
  report[['accuracy']] <- accuracy
  report[['kappa']] <- (accuracy - expAccuracy) / (1 - expAccuracy)
  report[['AUPRC_w']] <- mean(areas)
  report[['AUPRC_wc']] <- weighted.mean(areas, w)
  report[['multiLogLoss']] <- multiLogLoss
  report[['cm']] <- cm
  report[['confusionMatrix']] <- confusionMatrix
  report[['nonNA_PR']] <- nonNA_PR
  report[['PR_ROC']] <- PR_ROC
  
  return(report)
  
}




cal_class_ROCs<-function
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
  
  FPRfunc<-function(df){
    ans<-vector();
    for(i in 1:nrow(df)){
      ans<-append(ans, df[i,"FP"]/(df[i,"TN"]+df[i,"FP"]));
    }
    ans;
  }
  
  TPRfunc<-function(df){
    ans<-vector();
    for(i in 1:nrow(df)){
      ans<-append(ans, df[i,"TP"]/(df[i,"TP"]+df[i,"FN"]));
    }
    ans;
  }
  
  FPR<-FPRfunc(prsAll)
  TPR<-TPRfunc(prsAll)
  prsAll2<-cbind(prsAll, data.frame(TPR=TPR, FPR=FPR));
}

cal_class_PRs<-function
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
}

plot_class_ROCs<-function
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

  FPRfunc<-function(df){
    ans<-vector();
    for(i in 1:nrow(df)){
      ans<-append(ans, df[i,"FP"]/(df[i,"TN"]+df[i,"FP"]));
    }
    ans;
  }
  
  TPRfunc<-function(df){
    ans<-vector();
    for(i in 1:nrow(df)){
      ans<-append(ans, df[i,"TP"]/(df[i,"TP"]+df[i,"FN"]));
    }
    ans;
  }
  
  FPR<-FPRfunc(prsAll)
  TPR<-TPRfunc(prsAll)
  prsAll2<-cbind(prsAll, data.frame(TPR=TPR, FPR=FPR));

  ggplot(data=prsAll2, aes(x=as.numeric(as.vector(FPR)), y=as.numeric(as.vector(TPR)))) + geom_point(size = .5, alpha=.5) +  geom_path(size=.5, alpha=.75) +
  theme_bw() + xlab("False Positive Rate") + ylab("True Positive Rate") + facet_wrap( ~ ctype, ncol=4) +
  theme(axis.text = element_text(size=5)) + ggtitle("Classification performance")
}



#' @export
plot_metrics <- function(assessed){

 metric <- matrix(0, ncol = 2, nrow = 1)
 colnames(metric) <- c("cohen's kappa", "mean_AUPRC")
 rownames(metric) <- "value"
 metric[,1:2] <- c(assessed$kappa,assessed$AUPRC_w)

 df = as.data.frame(t(metric))
 df$eval = rownames(df)
 
 ggplot(df, aes(x = eval, y = value, fill = eval)) + geom_bar(stat = "identity") + 
  xlab("") + ylab("") + ylim(0, 1) + scale_fill_brewer(palette = "Set2") + 
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 8)) + 
  theme_bw() + theme(legend.position = "none")


}

#' @export
plot_PRs <- function(assessed, collapse=F){
  if(collapse){
    ggplot(data = assessed$nonNA_PR, aes(x = as.numeric(as.vector(recall)),y = as.numeric(as.vector(precision)),colour =ctype)) + 
      geom_point(size = 0.5, alpha = 0.5) + 
      geom_path(size = 0.5, alpha = 0.75) + 
      theme_bw() + xlab("Recall") + ylab("Precision") + 
      theme(axis.text = element_text(size = 5)) +
      ggtitle("Classification performance_PR Curve")
  }else{
    ggplot(data = assessed$nonNA_PR, aes(x = as.numeric(as.vector(recall)), y = as.numeric(as.vector(precision)))) + 
      geom_point(size = 0.5,alpha = 0.5) + 
      geom_path(size = 0.5, alpha = 0.75) + 
      theme_bw() + xlab("Recall") + ylab("Precision") + 
      facet_wrap(~ctype, ncol = 4) + 
      theme(axis.text = element_text(size = 5)) + 
      ggtitle("Classification performance_PR Curve")
  }
}

#' report sensitivity and precision of a scn score for a given cell type
#' @param score   a scn score for a given cell, numerical value
#' @param celltype  the cell type assigned by SCN
#' @param matrix  the helout assessment matrix obtained in SCN
#' @return a list of queried score, cell type, and PR
#' @export 
scn_calibration <- function(score, celltype, matrix){

  #input qc
  if(score > 1 || score < 0){
    print("Score is not within range! Please re-enter the score.")
  }
  
  if(!celltype %in% unique(matrix$ctype)){
    print("Cell type is not in the assessment matrix. Please re-enter celltype or check to see if there are cells left in heldout data to assess this cell type of interest.")
  }
  
  #figure out the lower bound and upper bounds
  score_rounded = round(score, digits = 2)
  if(score_rounded > score){
    score_lowerb = score_rounded - 0.005
    score_upperb = score_rounded
  }else{
    score_lowerb = score_rounded
    score_upperb = score_rounded + 0.005
  }
  confInt = c(score_lowerb, score_upperb)
  precision = c()
  recall = c()
  
  #access 
  precision = round(matrix[matrix$thresh %in% confInt & matrix$ctype == celltype, "precision"], digits = 3)
  recall = round(matrix[matrix$thresh %in% confInt & matrix$ctype == celltype, "recall"], digits = 3)
  
  print(paste("SCN score of", score, "for cell type", celltype, "has precision of", precision[1], "~", precision[2], 
              "and sensitivity of", recall[1], "~", recall[2], sep = " "))
  
  return(list(score = score, celltype=celltype, precision = precision, recall = recall))
  
}

plot_multiAssess <- function(assessed, method = "tsp_rf", ylimForMultiLogLoss = x){

 metric <- matrix(0, ncol = 5, nrow = 1)
  colnames(metric) <- c("cohen's kappa", "accuracy", "multiLogLoss", "mean_AUPRC","weighted-AUPRC")
  rownames(metric) <- "value"
  metric[,1:5] <- c(assessed$kappa, assessed$accuracy, assessed$multiLogLoss, assessed$AUPRC_w, assessed$AUPRC_wc) 
  metric <- as.data.frame(metric)

  p1<-ggplot(metric, aes(x="cohen's kappa", y = metric[1,1])) + geom_bar(stat="identity") +xlab("") + ylab("") + theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +  ylim(0,1) + theme(legend.position="none")

  p2<-ggplot(metric, aes(x="accuracy", y = metric[1,2])) + geom_bar(stat="identity") +xlab("") + ylab("") + theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) + ylim(0,1) + theme(legend.position="none")

  p3<-ggplot(metric, aes(x="multiLogLoss", y = metric[1,3])) + geom_bar(stat="identity") +xlab("") + ylab("") + theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) + ylim(0,ylimForMultiLogLoss)+ theme(legend.position="none")

  p4<-ggplot(metric, aes(x="mean_AUPRC", y = metric[1,4])) + geom_bar(stat="identity") +xlab("") + ylab("") + theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) + ylim(0,1) + theme(legend.position="none")

  p5<-ggplot(metric, aes(x="weight_AUPRC", y = metric[1,5])) + geom_bar(stat="identity") +xlab("") + ylab("") + theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) + ylim(0,1) + theme(legend.position="none")


  p6 <- ggplot(data=assessed$nonNA_PR, aes(x=as.numeric(as.vector(recall)), y=as.numeric(as.vector(precision)))) + geom_point(size = .5, alpha=.5) +  geom_path(size=.5, alpha=.75) +
    theme_bw() + xlab("Recall") + ylab("Precision") + facet_wrap( ~ ctype, ncol=4) +
    theme(axis.text = element_text(size=5)) + ggtitle("Classification performance_PR Curve")

  (p1 | p2 | p4 | p5 | p3) /
      p6
  
}

#select balanced training sample
selectTrain <- function(stDat, nCells, dLevel, sample_name){
  classes <- unique(stDat[,dLevel])
  goodClasses <- vector()
  newstTrain <- data.frame()
  
  for (i in 1: length(classes)){
    if(nrow(stDat[which(stDat[,dLevel] == classes[i]),]) > nCells){
      goodClasses <- c(goodClasses, classes[i])
      rowindex <- sample(rownames(stDat[which(stDat[,dLevel] == classes[i]),]), nCells)
      newstTrain <- rbind(newstTrain, stDat[rowindex,])
    }
  }
  
  return(newstTrain)
} 
