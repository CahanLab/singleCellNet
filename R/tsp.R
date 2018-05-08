filterTrain <- function(stTrain, dLevel = "description1", minCells = 30){
  goodGrps<-names(which(table(stDat_sub[,dLevel])>=minCells))
  newTrain<-data.frame()
  for(ggood in goodGrps){
    newTrain<-rbind(newTrain, stDat_sub[stDat_sub[,dLevel]==ggood,])
  }
  return(newTrain)
}

makeMultiTSP<-function(expDat, sampTab, topN=NULL){
  cell.types<-unique(as.vector(sampTab$description1))

  classList<-list()
  for(cell.type in cell.types){
    cat(cell.type,"... ")
    t1<-tsp_prep(expDat, sampTab, cell.type)
    myscores<-convertToDF(t1)
    classList[[cell.type]]<-makeClassifier(myscores, topN=topN)
    cat("\n")
  }
  classList
}

# make a matrix indicating gene vs gene props
sub_tsp_prep<-function(expDat, stX){
	#stX<-stTrain[stTrain$description1==desc,]
	expX<-expDat[,rownames(stX)]
	ngenes<-nrow(expX)
	ans<-matrix(0, nrow=ngenes, ncol=ngenes)
	for(i in 1:ngenes){
		for(j in 1:ngenes){
			if(j>i){
				ans[i,j]<- length( which(expX[i,]>expX[j,]) )
			}
		}
	}
	rownames(ans)<-rownames(expX)
	colnames(ans)<-rownames(expX)
	ans<-ans/ncol(expX)
}

tsp_prep<-function(expDat, stX, desc){

  st1<-stX[stX$description1==desc,]
  a1<-sub_tsp_prep(expDat, st1)

  st2<-stX[stX$description!=desc,]
  a2<-sub_tsp_prep(expDat, st2)

  # want to maintain directionality for later classification
  a1-a2
  
}

convertToDF<-function( tspRes ){
  n<-nrow(tspRes)
  genes<-rownames(tspRes)

  gene1<-rep("", length=(n**2)/2)
  gene2<-rep("", length=(n**2)/2)
  score<-rep(0, length=(n**2)/2)

  ci<-0
  for(i in 1:n){
    for(j in 1:n){
      if(j>i){
        ci<-ci+1

        gene1[ci]<-genes[i]
        gene2[ci]<-genes[j]
        score[ci]<-tspRes[i,j]
      }
    }
  }
  data.frame(gene1=gene1[1:ci], gene2=gene2[1:ci], score=score[1:ci], stringsAsFactors=FALSE)
}

# returns a list of gene pairs, order indicates class
makeClassifier<-function(tspScores, topN=NULL){

  if(is.null(topN)){
    xscores<-abs(tspScores$score)
    maxVal<-max(xscores)
    xi<-which(xscores==maxVal)
    tx<-tspScores[xi,]
  }
  else{
    tx<-tspScores[order( abs(tspScores$score), decreasing=TRUE),][1:topN,]
  }

  cat("Number of TSP: ",nrow(tx))
  ans<-list()
  for(i in 1:nrow(tx)){
    if(tx[i,]$score>0){
      ans[[i]]<-c(tx[i,"gene1"], tx[i,"gene2"])
    }
    else{
      ans[[i]]<-c(tx[i,"gene2"], tx[i,"gene1"])
    }
  }
  ans
}
#predict classs
predictClasses<-function(cList, expDat){
  ans<-matrix(0, nrow=length(cList), ncol=ncol(expDat))
  for(i in seq(length(cList))){
    ans[i,]<-predictClass(cList[[i]], expDat)
  }
  colnames(ans)<-colnames(expDat)
  rownames(ans)<-names(cList)
  ans
}

predictClass<-function(tspList, expDat){

   afunc<-function(aPair, expDat){
     expDat[aPair[1],]>expDat[aPair[2],]
   }

 
   tmpAns<-matrix(0, nrow=length(tspList),ncol=ncol(expDat))
   for(i in seq(length(tspList))){
     tmpAns[i,]<- as.numeric(afunc(tspList[[i]],expDat=expDat))
   }
   ssums<-apply(tmpAns, 2, mean)
   ssums
 }

#the following 3 functions can be found in CellNet github page
cn_classAssess<-function
(ct_scores,
 stVal, 
 classLevels="description1",
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