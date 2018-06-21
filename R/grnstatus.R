# Patrick Cahan (C) 2017
# patrick.cahan@gmail.com


#' barplot this specific GRN
#'
#' barplot this specific GRN
#' @param qScores queryScores
#' @param control scores
#' @param stQuery sample table
#' @param dLevelQ dLevel of query samples
#' @param snName name of subnet to plot establishment level 
#' @param ctrSamps names of samples in training data
#' @param bOrder order of bars
#' @param sidCol sample id colname
#'
#' @return gpplot barplot
cn_barplot_grnSing_base<-function### 
(qScores, 
 ctrlScores,
 stQuery,
 dLevelQ,
 snName,
 ctrSamps,
 bOrder,
 sidCol="sample_id"
 ){
    
  # convert into a data.frame
  aa<-scn_extract_SN_DF(qScores, stQuery,dLevelQ, rnames=snName, sidCol=sidCol);
  aa3<-scn_reduceMatLarge(aa, "score", "description", "subNet");
  aa3<-cbind(aa3, src=rep('query', nrow(aa3)));
  tmpAns<-data.frame();
  for(ctrSamp in ctrSamps){
    xxx<-ctrlScores[ctrlScores$grp_name==ctrSamp & ctrlScores$subNet==snName,];
    xxx$grp_name<-paste(xxx$grp_name, "_train", sep='');
    tmpAns<-rbind(tmpAns, xxx);
  }
  tmpAns<-cbind(tmpAns, src=rep("train", nrow(tmpAns)));
  
  if(is.null(bOrder)){
    bOrder<-c(as.vector(tmpAns$grp_name), as.vector(aa3$grp_name));
    ##bOrder<-aa3$grp_name[order(aa3$mean, decreasing=TRUE)];    
  }
  
  aa3<-rbind(aa3, tmpAns);
  aa3$grp_name<-factor(aa3$grp_name, bOrder);
  
  ##
  # convert is.na(stdev) -> 0
  xi<-which(is.na(aa3$stdev));
  if(length(xi)>0){
    aa3[xi,'stdev']<-0;
  }
  # # 
  ans<-  ggplot(na.omit(aa3), aes(x=grp_name, y=mean, fill=src)) +
    geom_bar(width=.75,position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=mean-stdev, ymax=mean+stdev),width=.2,position=position_dodge()) +
    scale_fill_brewer(palette = "Paired")  +
    theme_bw() +
    theme(text = element_text(size=8), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
    ggtitle(snName) + theme(axis.title.x = element_blank())+ ylab("GRN status")
  ans;
}

zscore<-function
### compute zscore
(x,
 ### numeric vector
 meanVal, 
 ### mean of distribution to compute zscore of x against 
 sdVal
 ### standard deviation of distribution to compute zscore of x agains
 ){ 
  (x-meanVal)/sdVal;
  ### zscore
}

#' Normalize grn status as compared to training data
#'
#' Divide the query scores by the mean values in the training data.
#' @param ctrlScores a list of subnet->mean value, all subnets
#' @param queryScores a matrix, rownames = subnet names, all subnets
#' @param subNets a vector of subnets names to use
#'
#' @return normalized grn status matrix
#'
scn_normalizeScores<-function
(ctrlScores, 
 queryScores, 
 subNets 
){
  
  ans<-matrix(0, nrow=length(subNets), ncol=ncol(queryScores));
  rownames(ans)<-subNets
  #subNets<-rownames(queryScores);
  for(subNet in subNets){
    ### cat(subNet,"\n")
    ans[subNet,]<- queryScores[subNet,] / ctrlScores[[subNet]];
  }
  colnames(ans)<-colnames(queryScores);
  ans;
}


#' computes the raw score for a gene as xmax-abs(zscore).
#'
#' better values are higher
#' @param vect a vector of gene expression values for multiple samples
#' @param mmean mean value in training data
#' @param ssd standard deviation in training data
#'
#' @return transformed (but not normalized) GRN score
#'
scn_rawScore<-function
(vect,
 mmean,
 ssd,
 xmax=1e3
){
  zcs<-zscore(vect, mmean, ssd);
  ### xmax<-1000; # arbitrary, and corrected for later, but want some high enough that it should not be exceeded 
  if(length(which(abs(zcs)>xmax)>0)) { cat(".") }
  xmax-abs(zcs);

}


utils_reduceMat<-function
### reduce the data.matrix values by averaging and getting st dvs
(datFrame,
 ### data frame to reduce
 valCol,
 ### value column
 cName='ann'
 ### cname to reduce on e.g. ann
){
    
  mids<-unique(as.vector(datFrame[,cName]));
  means<-rep(0, length(mids));
  sds<-rep(1, length(mids));
  indexI<-rep(0, length(mids)); # index to extract other columns from original data.frame
  
  for(i in seq(length(mids))){
    mid<-mids[i]
    xi<-which(datFrame[,cName]==mid);
    tmpDat<-datFrame[xi,];
    means[i]<-mean(tmpDat[,valCol]);
    sds[i]<-sd(tmpDat[,valCol]);
  }  
  data.frame(grp_name=mids, mean=means, stdev=sds);
  ### df of grp_name, mean, sd
}


scn_reduceMatLarge<-function
### reduce a data.frame
(datFrame,
 ### df
 valCol="score",
 ### value to merge
 cName="description",
 ### column to reduce on 
 iterOver="subNet"
 ### iterate over
 ){
  
  iterOvers<-unique(as.vector(datFrame[,iterOver]));
  ans<-data.frame();
  for(io in iterOvers){
    #  cat(io,"\n");
    xi<-which(datFrame[,iterOver]==io);
    dfTmp<-datFrame[xi,];
    x<- utils_reduceMat(dfTmp,valCol=valCol,cName=cName);
    x<-cbind(x, subNet=rep(io, nrow(x)));
    ans<-rbind(ans, x);
  }
  ans;
  ### ans
}


#' returns a DF of: sample_id, description, ctt, subnet_name, score
#'
#' returns a DF of: sample_id, description, ctt, subnet_name, score
#' @param scores a matrix of subNet scores
#' @param sampTab sample table
#' @param dLevel column name of sampTab to group values by
#' @param rnames rownames to extract
#' @param sidCol sample identifier column name
#'
#' @return returns a DF of: sample_id, description, ctt, subnet_name, score
scn_extract_SN_DF<-function
(scores,
 sampTab,
 dLevel,
 rnames=NULL,
 sidCol="sample_id"
){
  
  if(is.null(rnames)){
    rnames<-rownames(scores);
    #cat("GOT NULL\n");
  }
  tss<-scores[rnames,];
  if(length(rnames)==1){
    tss<-t(as.matrix(scores[rnames,]));
    rownames(tss)<-rnames;
  #  cat(dim(tss),"\n")
  }
  nSamples<-ncol(tss);
  stTmp<-sampTab[colnames(tss),]; ####
  snNames<-rownames(tss);
  num_subnets<-length(snNames);    
  snNames<-unlist(lapply(snNames, rep, times=nSamples));
  sample_ids<-rep(as.vector(stTmp[,sidCol]), num_subnets);
  descriptions<-rep(as.vector(stTmp[,dLevel]), num_subnets);
    # myCtts<-rep(ctt, length(snNames));
  scores<-as.vector(t(tss));
  data.frame(sample_id=sample_ids,
             description=descriptions,
             #         ctt=myCtts,
             subNet = snNames, 
             score=scores);
  ### data.frame
}


#' GRN status
#'
#' Calculates the status of all GRNs in query samples as compared to training data for
#' @param expDat query expression matrix
#' @param subList of ct => genes
#' @param tVals tvals
#' @param classList classList
#' @param minVals minVals
#' @param classWeight class weight
#' @param exprWeight  expression weight
#' @return grn scores (not normalized)
scn_netScores<-function
### return the GRN establishment score for a given expression matrix
(expDat, 
 genes, 
 tVals, 
 ctt, 
 classList=NULL, 
 classWeight=FALSE, 
 exprWeight=TRUE,
 xmax=1e3
){
  cat(ctt,"\n")
  aMat<-matrix(0, nrow=length(genes), ncol=ncol(expDat));
  rownames(aMat)<-genes;
  
  weights<-rep(1, length(genes));
  names(weights)<-genes;
  
  #otherCTs<-setdiff(names(tVals), ct)
  
  cat(dim(aMat),"\n")
  if(exprWeight){
    meanVect<-unlist(tVals[[ctt]][['mean']][genes]);
    weights<-(2**meanVect)/sum(2**meanVect);
    if(FALSE){
      bDifs<-minDif(tVals, genes, ctt)
    # set to zero neg vals
      bDifs[bDifs<0]<-0
      weights<-bDifs/sum(bDifs)
    }
  }
    
  if(classWeight){
    classImp<-classList[[ctt]]$importance[genes,1];
    ### 04-19-17
    ###classImp<-classImp/sum(classImp)
    weights<-weights*classImp;
  }
  
  for(gene in genes){
   ### cat("***",gene,"\n")
    ###zzs<-as.matrix(cn_rawScore(expDat[gene,], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]])[1,])


    zzs<-scn_rawScore(expDat[gene,], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]], xmax=xmax)
    aMat[gene,]<-zzs;

  }
  xscores<-apply(aMat, 2, weighted.mean, w=weights);
  xscores;
}

#' GRN status
#'
#' Calculates the GRN status in query samples as compared to training data
#' @param expDat query expression matrix
#' @param subList of ct => genes
#' @param tVals list of ctt->list(means->genes, sds->genes)
#' @param classList class list
#' @param minVals  min vals
#' @param classWeight classweght
#' @param exprWeight  expression weight
#' @return GRN scores 
#'
scn_score<-function
(expDat,
 subList, 
 tVals,
 classList=NULL,
 minVals=NULL, 
 classWeight=FALSE,
 exprWeight=TRUE,
 xmax=1e3
){
  #nSubnets<-sum(sapply(subList, length));
  nSubnets<-length(subList);
  ans<-matrix(0, nrow=nSubnets, ncol=ncol(expDat));
  ctts<-names(subList);
  rnames<-vector();
  rIndex<-1;
  for(ctt in ctts){
     cat(ctt,"\n");
    genes<-subList[[ctt]];
    # 06-06-16 -- added to allow for use of GRNs defined elsewhere
    genes<-intersect(genes, rownames(expDat));
    #    snNames<-names(subnets);
    #    rnames<-append(rnames, snNames);    
    #    for(sName in snNames){
    ans[rIndex,]<-scn_netScores(expDat, genes, tVals=tVals, ctt=ctt,classList, classWeight=classWeight,exprWeight=exprWeight, xmax=xmax);
    rnames<-append(rnames, ctt);
    rIndex<-rIndex+1;
    #   }
  }
  rownames(ans)<-rnames;
  colnames(ans)<-colnames(expDat);
  if(!is.null(minVals)){
    minVals<-minVals[rownames(ans)];
    ans<-ans-minVals;
  }
  ans;
}




#
# functions to enable GRN status metric
#

#' Figure out normalization factors for GRNs, and norm training data
#'
#' Exactly that.
#' @param expTrain expression matrix
#' @param stTrain  sample table
#' @param subNets named list of genes, one list per CTT, tct=>gene vector
#' @param classList list of classifiers
#' @param dLevel column name to group on
#' @param tVals seful when debugging
#' @param classWeight weight GRN status by importance of gene to classifier
#' @param exprWeight weight GRN status by expression level of gene?
#' @param sidCol sample id colname
#' 
#' @return list of trainingScores, normVals, raw_scores, minVals, tVals=tVals
#' @export
scn_trainNorm<-function # 
(expTrain,
 stTrain,
 subNets,
 classList = NULL, 
 dLevel = "description1",
 tVals=NULL,
 classWeight=FALSE,
 exprWeight=FALSE,
 sidCol='sample_id',
 xmax=1e3,
 predSD=FALSE
){

  if(is.null(tVals)){
    tVals<-scn_make_tVals(expTrain, stTrain, dLevel)
  }

  ctts<-as.vector(unique(stTrain[,dLevel]));
  scoreList<-list();
  normList<-list(); # a list of ctt->subnet->mean value
  minVect<-vector(); # a list of ctt->subnet->min value, used to shift raw grn est scores
  
  cat("calculating GRN scores on training data ...\n");
  tmpScores<-scn_score(expTrain, subNets, tVals, classList, minVals=NULL, classWeight=classWeight, exprWeight=exprWeight, xmax=xmax)


  minVect<-apply(tmpScores, 1, min);


  names(minVect)<-rownames(tmpScores);
  
  # shift the raw scores so that min=0;

  tmpScores<-tmpScores - minVect;
  cat("norm factors\n");
  for(ctt in ctts){
    # determine nomalization factors
    ##snets<-names(subNets[[ctt]]);
    snets<-ctt;

    scoreDF<-scn_extract_SN_DF(tmpScores, stTrain, dLevel, snets, sidCol=sidCol);
    scoreDF<-scn_reduceMatLarge(scoreDF, "score", "description", "subNet");
    xdf<-scoreDF[which(scoreDF$grp_name==ctt),];
    tmpSNS<-as.list(xdf$mean);
    names(tmpSNS)<-xdf$subNet;
    normList[names(tmpSNS)]<-tmpSNS;      
  }
  
  # normalize training scores
  nScores<-scn_normalizeScores(normList, tmpScores, rownames(tmpScores));

  scoreDF<-scn_extract_SN_DF(nScores, stTrain, dLevel, sidCol=sidCol);

  scoreDF<-scn_reduceMatLarge(scoreDF, "score", "description", "subNet");

  list(trainingScores=scoreDF,
       normVals=normList,
       raw_scores=tmpScores,
       minVals=minVect,
       tVals=tVals);
}

#' Estimate gene expression dist in CTs
#'
#' Calculate mean and SD 
#' @param expDat training data 
#' @param sampTab, ### training sample table 
#' @param dLevel="description1", ### column to define CTs
#' @param predictSD=FALSE ### whether to predict SD based on expression level
#'
#' @return tVals list of ct->mean->named vector of average gene expression, ->sd->named vector of gene standard deviation
scn_make_tVals<-function### estimate gene expression dist in CTs
(expDat, ### training data 
 sampTab, ### training sample table 
 dLevel="description1" ### column to define CTs
){
  
    # Note: returns a list of dName->gene->mean, sd, where 'dName' is a ctt or lineage 
    # make sure everything is lined up
    expDat<-expDat[,rownames(sampTab)];
    tVals<-list();
    dNames<-unique(as.vector(sampTab[,dLevel]));
    allGenes<-rownames(expDat);
    for(dName in dNames){
      #cat(dName,"\n");
      xx<-which(sampTab[,dLevel]==dName);
      sids<-rownames(sampTab[xx,]);
      xDat<-expDat[,sids];
      means<-apply(xDat, 1, mean);
      sds<-apply(xDat, 1, sd);
      tVals[[dName]][['mean']]<-as.list(means);
      tVals[[dName]][['sd']]<-as.list(sds);
    }
    ans<-tVals;
  ans;
}
