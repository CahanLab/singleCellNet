############################################################Calculating C/T-specific GRN status

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
netScores<-function (expDat, genes, tVals, ctt, classList=NULL, classWeight=TRUE, classWeightVal = 3, exprWeight=TRUE, exprWeightVal = 3, xmax=1e3){
  cat(ctt,"\n")
  aMat<-matrix(0, nrow=length(genes), ncol=ncol(expDat));
  rownames(aMat)<-names(genes);
  
  weights<-rep(1, length(genes));
  names(weights)<-names(genes);
  
  #otherCTs<-setdiff(names(tVals), ct)
  
  cat(dim(aMat),"\n")
  if(exprWeight){
    meanVect<-unlist(tVals[[ctt]][['mean']][names(genes)]);
    weights<-(exprWeightVal*meanVect)/sum(exprWeightVal*meanVect); # also arbritary value on the weight you are putting on the expression
  }
  
  if(classWeight){ #TODO modify this to fit the gene pairs
    #classImp<-classList[[ctt]]$importance[genes,1];
    
    classImp = weights
    for(gene in names(classList[[ctt]])) {
      if(gene %in% names(classImp)) {
        classImp[gene] = classList[[ctt]][gene] + classWeightVal # arbritary value
      }
    }
    
    ### 04-19-17
    ###classImp<-classImp/sum(classImp)
    weights<-weights*classImp;
  }
  
  for(gene in names(genes)){
    
    ### cat("***",gene,"\n")
    ###zzs<-as.matrix(cn_rawScore(expDat[gene,], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]])[1,])
    
    zzs<-rawScore(expDat[gene,], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]], xmax=xmax, reg_type = genes[gene])
    
    aMat[gene,] = zzs;
  }
  
  print(dim(aMat))
  xscores<-apply(aMat, 2, weighted.mean, w=weights);
  xscores;
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
make_tVals<-function (expDat, sampTab, dLevel="description1", predictSD=FALSE){
  
  if(predictSD){
    ans<-make_tVals_predict(expDat, sampTab, dLevel);
  }
  else{
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
  }
  ans;
}

#' make_tVals_predict
#'
#' predicts SD based on mean expression using linear regression
#' @param expDat training data
#' @param sampTab training sample table
#' @param dLevel="description1" column to define CT
#'
#' @return tVals list of ct->mean->named vector of average gene expression, ->sd->named vector of gene standard deviation
make_tVals_predict<-function(expDat, sampTab, dLevel="description1"){
  # Note: returns a list of dName->gene->mean, sd, where 'dName' is a ctt or lineage
  # make sure everything is lined up
  expDat<-expDat[,rownames(sampTab)];
  tVals<-list();
  dNames<-unique(as.vector(sampTab[,dLevel]));
  allGenes<-rownames(expDat);
  
  # make a model to predict SD given average expression level across all samples
  sdT<-apply(expDat, 1, sd);
  mT<-apply(expDat, 1, mean);
  myModel<-lm(sdT~mT);
  for(dName in dNames){
    xx<-which(sampTab[,dLevel]==dName);
    sids<-rownames(sampTab[xx,]);
    xDat<-expDat[,sids];
    means<-apply(xDat, 1, mean);
    sds<-predict(myModel, data.frame(mT=means));
    tVals[[dName]][['mean']]<-as.list(means);
    tVals[[dName]][['sd']]<-as.list(sds);
  }
  tVals;
}

#' Min Difference
#'
#' computes mean gene A in CT 1 - mean gene A in CT 2, where CT2 has the non CT1 max value. does this for genes
#' @param tVals tVals
#' @param genes vector of gene names
#' @param ct ct to compare to
#' @return vector of differences
minDif<-function(tVals, genes, ct){
  octs<-setdiff(names(tVals), ct)
  qq<-lapply(tVals[octs], "[[", "mean")
  ##tVals[[ct]][["mean"]][[gene]]#-max(unlist(lapply(qq, "[[", gene)))
  tmpMat<-matrix(unlist(lapply(qq, "[", genes)), nrow=length(genes))
  rownames(tmpMat)<-genes
  maxes<-apply(tmpMat, 1, max)
  unlist(tVals[[ct]][["mean"]][genes])-maxes
}

#' computes the raw z score for a gene as xmax-abs(zscore).
#'
#' better values are higher.
#' @param vect a vector of gene expression values for multiple samples
#' @param mmean mean value in training data
#' @param ssd standard deviation in training data
#' @param reg_type the type of regulation (up or down) for the specific subnetwork
#' @return transformed (but not normalized) GRN score
#'
rawScore<-function(vect, mmean, ssd, xmax=1e3, reg_type){
  zcs = zscore(vect, mmean, ssd);
  
  if(as.numeric(reg_type) == 1) { # if the
    #zcs[zcs > 0 & zcs < 1] = 0 # change this to revert back to original place
    zcs[zcs > 0] = 0
    return(xmax - abs(zcs))
  }
  else {
    #zcs[zcs < 0 & zcs > -1] = 0 # change this to revert back to original place
    zcs[zcs < 0] = 0
    return(xmax - abs(zcs))
  }
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
#' @export
grnScore<-function(expDat, subList, tVals, classList=NULL, minVals=NULL, classWeight=FALSE, exprWeight=TRUE, xmax=1e3){
  #nSubnets<-sum(sapply(subList, length));
  if(class(expDat) != "matrix") {
    expDat = as.matrix(expDat)
  }
  nSubnets<-length(subList);
  ans<-matrix(0, nrow=nSubnets, ncol=ncol(expDat));
  ctts<-names(subList);
  rnames<-vector();
  rIndex<-1;
  for(ctt in ctts){
    cat(ctt,"\n");
    genes<-subList[[ctt]];
    # 06-06-16 -- added to allow for use of GRNs defined elsewhere
    genes<-genes[intersect(names(genes), rownames(expDat))]; # only select the genes that are
    #    snNames<-names(subnets);
    #    rnames<-append(rnames, snNames);
    #    for(sName in snNames){
    ans[rIndex,]<-netScores(expDat, genes, tVals=tVals, ctt=ctt,classList=classList, classWeight=classWeight,exprWeight=exprWeight, xmax=xmax);
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

#' Normalize grn status as compared to training data
#'
#' Divide the query scores by the mean values in the training data.
#' @param ctrlScores a list of subnet->mean value, all subnets
#' @param queryScores a matrix, rownames = subnet names, all subnets
#' @param subNets a vector of subnets names to use
#'
#' @return normalized grn status matrix
#' @export
normalizeScores<-function(ctrlScores, queryScores, subNets){
  
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
#' @param xmax the maximum raw score that a sample could receive per gene
#' @param meanNorm normalize raw scores based on the lowest mean in a category
#'
#' @return list of trainingScores, normVals, raw_scores, minVals, tVals=tVals
#' @export
trainNorm<-function (expTrain, stTrain, subNets, classList = NULL,  dLevel = "description1", 
                     tVals=NULL, classWeight=TRUE, exprWeight=FALSE, sidCol='sample_id', xmax=1e3, meanNorm = FALSE){
  
  if(is.null(tVals)){
    tVals<-make_tVals(expTrain, stTrain, dLevel)
  }
  
  ctts<-as.vector(unique(stTrain[,dLevel]));
  scoreList<-list();
  normList<-list(); # a list of ctt->subnet->mean value
  minVect<-vector(); # a list of ctt->subnet->min value, used to shift raw grn est scores
  
  cat("calculating GRN scores on training data ...\n");
  tmpScores<-score(expTrain, subNets, tVals, classList, minVals=NULL, classWeight=classWeight, exprWeight=exprWeight, xmax=xmax)
  
  
  if(meanNorm == TRUE) {
    train_meanScores = meanTraining(tmpScores, stTrain, dLevel, sidCol)
    minVect<-apply(train_meanScores, 1, min);
    names(minVect)<-rownames(train_meanScores);
    
  } else {
    minVect<-apply(tmpScores, 1, min);
    names(minVect)<-rownames(tmpScores);
    
  }
  
  
  # shift the raw scores so that min=0;
  tmpScores<-tmpScores - minVect;
  cat("norm factors\n");
  for(ctt in ctts){
    # determine nomalization factors
    ##snets<-names(subNets[[ctt]]);
    snets<-ctt;
    
    scoreDF<-extract_SN_DF(tmpScores, stTrain, dLevel, snets, sidCol=sidCol);
    scoreDF<-reduceMatLarge(scoreDF, "score", "description", "subNet");
    xdf<-scoreDF[which(scoreDF$grp_name==ctt),];
    tmpSNS<-as.list(xdf$mean);
    names(tmpSNS)<-xdf$subNet;
    normList[names(tmpSNS)]<-tmpSNS;
  }
  
  # normalize training scores
  nScores<-normalizeScores(normList, tmpScores, rownames(tmpScores));
  
  scoreDF<-extract_SN_DF(nScores, stTrain, dLevel, sidCol=sidCol);
  
  scoreDF<-reduceMatLarge(scoreDF, "score", "description", "subNet");
  
  list(trainingScores=scoreDF,
       normVals=normList,
       raw_scores=tmpScores,
       minVals=minVect,
       tVals=tVals);
}

#' @title calculate the mean GRN scores across categories
#' @description calculate the mean GRN scores across categories
#' @param grnScores the GRN scores calculated through score
#' @param stTrain the sample table used for training
#' @param dLevel the name of the column with all the categories
#' @return an averaged GRN score matrix
meanTraining <- function(grnScores, stTrain, dLevel, sidCol) {
  rownames(stTrain) = as.vector(stTrain[, sidCol])
  # return matrix
  meanMatrix = matrix(0, nrow = nrow(grnScores), ncol = length(unique(stTrain[, dLevel])))
  
  rownames(meanMatrix) = rownames(grnScores)
  colnames(meanMatrix) = unique(stTrain[, dLevel])
  
  for(cancerName in unique(stTrain[, dLevel])) {
    tempStTrain = stTrain[stTrain[, dLevel] == cancerName, ]
    tempGRNscores = grnScores[, rownames(tempStTrain)]
    
    meanMatrix[, cancerName] = apply(tempGRNscores, 1, mean)
  }
  
  return(meanMatrix)
}
#' Z scores
#' Figure out Z score given mean and standard deviation
#' @param x query score
#' @param meanVal mean of the distribution
zscore<-function(x,meanVal,sdVal){
  
  (x-meanVal)/sdVal;
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
extract_SN_DF<-function(scores, sampTab, dLevel, rnames=NULL, sidCol="sample_id"){
  
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

#' reduce large matrix
#' reduce large matrix from extract_SN_DF
#'
#' @param datFrame the result from extract_SN_DF
#' @param valCol the column name of the score
#' @param cName the column name of the different groups
#' @param iterOver the column name of the subnetworks
#'
#' @return reduced large matrix
reduceMatLarge<-function (datFrame, valCol="score", cName="description", iterOver="subNet"){
  
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

#' @title  Reduce data matrix
#' @description reduce the data.matrix values by averaging and getting st dvs
#'
#' @param datFrame the dataframe
#' @param valCol the column with scores
#' @param cName column with groups
#'
#' @return df of grp_name, mean, sd
utils_reduceMat<-function(datFrame, valCol, cName='ann'){
  
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
}

#' @title Get Query GRN status
#' @description Get the GRN status of query samples
#'
#' @param expQuery logRanked query expression matrix
#' @param expTrain logRanked training expression matrix
#' @param stTrain sample table of training expression matrix
#' @param dLevel the name of the column with cancer types
#' @param sidCol the name of the column with sample IDs
#' @param grn_return the grn list that is returned from \code{\link{scn_makeGRN}}
#' @param trainNorm normalization statistics from \code{\link{trainNorm}}. If you are using pre-calculated normalization statistics please make sure all the parameters are the same for applying to query and calculating normalization
#' @param classifier_return the classifier_return list that is returned from \code{\link{broadClass_train}}
#' @param classWeight boolean indicating whether to take the importance of the classification into status calculation
#' @param exprWeight boolean indicating whether to take the weight of gene expression into status calculation
#' @param prune boolean indicating whether to select exclusive genes for processing classification gene importance
#' @param predSD a parameter for calculating normalization statistics from training data
#' @return a matrix indicating the GRN status
#' @export
queryGRNstatus <- function(expQuery, expTrain, stTrain, dLevel, sidCol, grn_return, trainNorm = NULL, classifier_return,  classWeight = TRUE, exprWeight = FALSE, prune = TRUE, xmax = 1e3, predSD=FALSE) {
  cnProc = classifier_return$cnProc
  
  trainNorm_prior = trainNorm
  geneImportance = processImportance(classifier = cnProc$classifier, xpairs = classifier_return$xpairs_list, prune = prune)
  
  if(is.null(trainNorm) == TRUE) {
    trainNorm = trainNorm(expTrain, stTrain, subNets=grn_return$ctGRNs$geneLists, classList = geneImportance, dLevel = dLevel, sidCol = sidCol, classWeight = classWeight, exprWeight = exprWeight)
    cat("Finished finding normalization statistics", "\n")
  }
  
  status_score = grnScore(expDat = expQuery,
                           subList = grn_return$ctGRNs$geneLists, tVals = trainNorm$tVals,
                           classList = geneImportance, minVals = trainNorm$minVals,
                           classWeight = classWeight, exprWeight = exprWeight,
                           xmax = xmax)
  normScoresQuery = normalizeScores(trainNorm$normVals, status_score, rownames(status_score))
  
  if(is.null(trainNorm_prior) == TRUE) {
    return(list("trainNorm" = trainNorm, "query_GRNstatus" = normScoresQuery))
  }
  else {
    return(normScoresQuery)
  }
}

#' @title score indvidual genes
#' @description score individual genes
#' @param queryMatrix query ranked expression matrix
#' @param ctt the name of the cancer type specific network
#' @param trainNorm the normalization statistics from \code{\link{trainNorm}}
#' @param grn_all the grn construction results from \code{\link{makeGRN}}
#' @param importantGenes a list of gene importances
#' @return a matrix with individual gene z score for each query sample, and the expression direction it should move to achieve ideal similarity with subnetwork. The order of the gene is ranked by importance from the classifier
#' @export
geneScores <- function(queryMatrix, ctt, trainNorm, grn_all, importantGenes) {
  tvals = trainNorm$tVals
  cgenes = grn_all$ctGRNs$geneLists[[ctt]]
  
  returnMatrix = matrix(nrow = length(cgenes), ncol = ncol(queryMatrix))
  rownames(returnMatrix) = names(cgenes)
  colnames(returnMatrix) = colnames(queryMatrix)
  
  for(gene in names(cgenes)) {
    returnMatrix[gene, ] = cgenes[[gene]] * (1000 - rawScore(queryMatrix[gene, ], tvals[[ctt]][["mean"]][[gene]], tvals[[ctt]][["sd"]][[gene]], reg_type = cgenes[[gene]], xmax = 1000))
    
  }
  
  genes_important = importantGenes[[ctt]]
  orderGenes = genes_important[intersect(names(genes_important), rownames(returnMatrix))]
  orderGenes = sort(orderGenes, decreasing = TRUE)
  
  notImportantGenes = rownames(returnMatrix)[!(rownames(returnMatrix) %in% names(orderGenes))]
  
  orderGenes = c(names(orderGenes), notImportantGenes)
  returnMatrix = returnMatrix[orderGenes, ]
  return(returnMatrix)
}





############################################################Calculating TF score

#' @title calculate TF scores
#' @description calculate the TF that one should target to yield target cell type
#'
#' @param expQuery the logRanked expression matrix for query data
#' @param subnetName the name of the target cell type
#' @param grnAll result of running \code{\link{scn_makeGRN}}
#' @param trainNorm result of running \code{\link{scn_trainNorm}}
#' @param classifier_return the classifier_return list that is returned from \code{\link{broadClass_train}}
#' @param classWeight whether to take the weight of the classifier into calculation
#' @param classWeightVal the value of the classifier weight
#' @param exprWeight whether to take the weight of the expression into calculation
#' @param exprWeightVal the value of the expression weight
#' @param correlationFactor the weight of correlation direction in the network
#' @param prune the pararmeter for pruning importance of gene pairs 
#' @param normTFscore boolean indicate whether to normalize TF scores
#'
#' @return matrix of TF scores and query samples
#' @export
tfScores <- function(expQuery, subnetName, grnAll, trainNorm, classifier_return, classWeight=TRUE, classWeightVal = 3, exprWeight=FALSE, exprWeightVal = 3, correlationFactor = 1, prune = TRUE, normTFscore = FALSE) {
  cnProc = classifier_return$cnProc
  
  classList = processImportance(classifier = cnProc$classifier, xpairs = classifier_return$xpairs_list, prune = prune)
  
  
  TF_targetList = grnAll[["ctGRNs"]][["tfTargets"]][[subnetName]]
  
  tfs = names(TF_targetList)
  netGenes = grnAll[["ctGRNs"]][["geneLists"]][[subnetName]]
  
  grnTable = grnAll$overallGRN$grnTable
  tVals = trainNorm$tVals
  
  # get the z scores matrix
  z_scoreMat = matrix(0, nrow = length(netGenes), ncol = ncol(expQuery))
  rownames(z_scoreMat) = names(netGenes)
  colnames(z_scoreMat) = colnames(expQuery)
  
  # get the tf scores matrix
  tf_scoreMat = matrix(0, nrow = length(tfs), ncol = ncol(expQuery))
  rownames(tf_scoreMat) = tfs
  colnames(tf_scoreMat) = colnames(expQuery)
  
  print(tfs)
  for(querySample in colnames(expQuery)) {
    xvals = as.vector(expQuery[names(netGenes), querySample])
    names(xvals) = names(netGenes)
    z_scoreMat[, querySample] = zscoreVect(genes = names(netGenes), xvals = xvals, tVals = tVals, ctt = subnetName)
    
  }
  
  # assign weights to each gene
  weights = rep(1, length(netGenes))
  names(weights) = names(netGenes)
  if(exprWeight){
    meanVect = unlist(tVals[[subnetName]][['mean']][names(netGenes)]);
    weights = (exprWeightVal*meanVect)/sum(exprWeightVal*meanVect); # also arbritary value on the weight you are putting on the expression
  }
  
  if(classWeight){ #TODO modify this to fit the gene pairs
    classImp = weights
    for(gene in names(classList[[subnetName]])) {
      if(gene %in% names(classImp)) {
        classImp[gene] = classList[[subnetName]][gene] + classWeightVal # arbritary value
      }
    }
    weights = weights*classImp;
  }
  
  
  for(sampleID in colnames(z_scoreMat)) {
    
    tfScore = vector() # the vector of tf z-scores
    # to calculate TFs
    cat("Calculating TF scores for", sampleID, "\n")
    
    for(i in seq(1, length(tfs))) {
      tf = tfs[i]
      targs = TF_targetList[[tf]]; #
      targs = intersect(targs, rownames(expQuery))
      
      temp_tfScore = calc_tfScores(tf, targs, sampleID, z_scoreMat, netGenes, weights, grnTable, correlationFactor)
      if(normTFscore == TRUE) {
        temp_tfScore = temp_tfScore - calc_normTfScore(tf, targs, sampleID, netGenes, weights, grnTable, correlationFactor)
        
      }
      tfScore = c(tfScore, temp_tfScore)
    }
    
    tf_scoreMat = tf_scoreMat[names(tfScore), ]
    tf_scoreMat[, sampleID] = tfScore
  }
  
  return(tf_scoreMat)
}

#' @title Compute mean Z scores of given genes
#' @description Compute the mean Z score of given genes in each sample
#'
#' @param genes the genes
#' @param xvals named vector
#' @param tVals mean and SD of genes in training data
#' @param ctt cell type
#'
#' @return a vector of Z scores
zscoreVect<-function (genes, xvals, tVals, ctt){
  ans<-vector();
  for(gene in genes){
    ans<-append(ans, zscore(xvals[gene], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]]));
  }
  ans;
}

#' Calculate the TF
#' Calculate the TF score given all the parameters
#'
#' @param tf TF
#' @param targs a list of targets for TF
#' @param sampleID the sample ID
#' @param z_scoreMat the z value matrix
#' @param netGenes genes selected for building the subnetwork
#' @param weights a vector of weights for each gene
#' @param grnTable the grn table with correlation values and signs
#' @param correlationFactor the weight of correlation between
calc_tfScores <- function(tf, targs, sampleID, z_scoreMat, netGenes, weights, grnTable, correlationFactor = 1) {
  
  if(z_scoreMat[tf, sampleID] > 0 & netGenes[tf] == 1) { # if the gene suppose to be higher and is above z
    temp_zscore = 1
  } else if(z_scoreMat[tf, sampleID] < 0 & netGenes[tf] == -1) { # if the gene is suppose to be lower and is below z
    temp_zscore = 1
  } else {
    temp_zscore = 1 + abs(z_scoreMat[tf, sampleID])
  }
  
  part1 = length(targs) * temp_zscore * weights[tf]
  
  part2 = 0
  
  # revamp this shit
  calculation_matrix = matrix(data = 1, nrow = length(targs), ncol = 8)
  rownames(calculation_matrix) = targs
  colnames(calculation_matrix) = c("TF_direction", "Target_direction","z_score", "z_mod", "weight", "corr", "correlation_factor", "total_score")
  
  calculation_matrix[, "TF_direction"] = netGenes[tf] # assign the TF direction
  calculation_matrix[, "Target_direction"] = netGenes[targs] # assign the target direction
  
  temp_zscore = z_scoreMat[targs, sampleID]
  calculation_matrix[, "z_score"] = as.vector(temp_zscore)
  
  # assign modification z score
  modIndex = (calculation_matrix[, "Target_direction"] == 1 & calculation_matrix[, "z_score"] < 0) # if z score is low and its suppose to be upregulated
  calculation_matrix[modIndex, "z_mod"] = abs(calculation_matrix[modIndex, "z_score"]) + 1
  
  modIndex = (calculation_matrix[, "Target_direction"] == -1 & calculation_matrix[, "z_score"] > 0) # if z score is high and its suppose to be downregulated
  calculation_matrix[modIndex, "z_mod"] = abs(calculation_matrix[modIndex, "z_score"]) + 1
  
  # assign weight
  calculation_matrix[, "weight"] = weights[targs]
  
  # assign correlation direction
  temp_grnTable = grnTable[grnTable$TF == tf, ]
  rownames(temp_grnTable) = as.vector(temp_grnTable$TG)
  temp_grnTable = temp_grnTable[targs, ]
  corr_sign = sign(temp_grnTable[, "corr"])
  calculation_matrix[, "corr"] = corr_sign
  
  # assign the correlation factor
  calculation_matrix[, "correlation_factor"] = -correlationFactor # default to be negative
  
  # if correlation between target gene and TF is positive
  # and if both target and TF are suppose to be upregulated or downregulated
  # it's do not penalize
  modIndex = (calculation_matrix[, "corr"] == 1 & calculation_matrix[, "Target_direction"] == calculation_matrix[, "TF_direction"]) # if z score is low and its suppose to be upregulated
  calculation_matrix[modIndex, "correlation_factor"] = abs(correlationFactor)
  
  # if correlation between target gene and TF is negative
  # and if both target and TF are in opposite direction
  # do not penalize
  modIndex = (calculation_matrix[, "corr"] == -1 & calculation_matrix[, "Target_direction"] != calculation_matrix[, "TF_direction"]) # if z score is low and its suppose to be upregulated
  calculation_matrix[modIndex, "correlation_factor"] = abs(correlationFactor)
  
  # calculate total score for individual gene
  calculation_matrix[, "total_score"] = calculation_matrix[, "correlation_factor"] * calculation_matrix[, "weight"] * calculation_matrix[, "z_mod"]
  
  part2 = sum(calculation_matrix[, "total_score"])
  
  TF_score = part1 + part2
  
  # if TF is suppose to be down regulated
  if(netGenes[tf] == -1) {
    return(-(TF_score))
  }
  else {
    return(TF_score)
    
  }
  
}

#' Calculate the TF normalization constant
#' Calculate the TF score given all the parameters
#'
#' @param tf TF
#' @param targs a list of targets for TF
#' @param sampleID the sample ID
#' @param netGenes genes selected for building the subnetwork
#' @param weights a vector of weights for each gene
#' @param grnTable the grn table with correlation values and signs
#' @param correlationFactor the weight of correlation between
calc_normTfScore <- function(tf, targs, sampleID, netGenes, weights, grnTable, correlationFactor = 1) {
  
  part1 = length(targs) * weights[tf]
  part2 = 0
  
  # revamp this shit
  calculation_matrix = matrix(data = 1, nrow = length(targs), ncol = 8)
  rownames(calculation_matrix) = targs
  colnames(calculation_matrix) = c("TF_direction", "Target_direction","z_score", "z_mod", "weight", "corr", "correlation_factor", "total_score")
  
  calculation_matrix[, "TF_direction"] = netGenes[tf] # assign the TF direction
  calculation_matrix[, "Target_direction"] = netGenes[targs] # assign the target direction
  
  # assign weight
  calculation_matrix[, "weight"] = weights[targs]
  
  # assign correlation direction
  temp_grnTable = grnTable[grnTable$TF == tf, ]
  rownames(temp_grnTable) = as.vector(temp_grnTable$TG)
  temp_grnTable = temp_grnTable[targs, ]
  corr_sign = sign(temp_grnTable[, "corr"])
  calculation_matrix[, "corr"] = corr_sign
  
  # assign the correlation factor
  calculation_matrix[, "correlation_factor"] = -correlationFactor # default to be negative
  
  # if correlation between target gene and TF is positive
  # and if both target and TF are suppose to be upregulated or downregulated
  # it's do not penalize
  modIndex = (calculation_matrix[, "corr"] == 1 & calculation_matrix[, "Target_direction"] == calculation_matrix[, "TF_direction"]) # if z score is low and its suppose to be upregulated
  calculation_matrix[modIndex, "correlation_factor"] = abs(correlationFactor)
  
  # if correlation between target gene and TF is negative
  # and if both target and TF are in opposite direction
  # do not penalize
  modIndex = (calculation_matrix[, "corr"] == -1 & calculation_matrix[, "Target_direction"] != calculation_matrix[, "TF_direction"]) # if z score is low and its suppose to be upregulated
  calculation_matrix[modIndex, "correlation_factor"] = abs(correlationFactor)
  
  # calculate total score for individual gene
  calculation_matrix[, "total_score"] = calculation_matrix[, "correlation_factor"] * calculation_matrix[, "weight"] * calculation_matrix[, "z_mod"]
  
  part2 = sum(calculation_matrix[, "total_score"])
  
  TF_score = part1 + part2
  
  # if TF is suppose to be down regulated
  if(netGenes[tf] == -1) {
    return(-(TF_score))
  }
  else {
    return(TF_score)
  }
  
}

#' @title compute context dependent zscores
#'
#' @description slightly modidied from JJ Faith et al 2007
#' @param corrMat correlation matrix
#'
#' @return matrix of clr zscores
mat_zscores<-function(corrMat){
  corrMat<-abs(corrMat);
  zscs_2<-round(scale(corrMat), 3);
  rm(corrMat);
  gc()
  zscs_2 + t(zscs_2);
}


# utility functions

#' @title logRank gene expression
#' @description rank the gene expression values and log10 the ranks
#' @param expTrain the expression matrix
#' @return the log10 rank of genes
#' @export
logRank <- function(expTrain, base = 2) {
  expTrain = apply(expTrain, FUN = rank, MARGIN = 2)
  
  if (base == 0) {
    return(expTrain)
  }
  else {
    return(log(expTrain, base = base))
    
  }
}

#' @title process classifier importance
#' @description process the importance scores of gene pairs into enriched genes
#' @param class_info the class_info object return from scn_train
#' @param prune blooln indicate whether to select genes exclusive to cancer category
#'
#' @return a list of genes with their importance divided by groups
#' @export
processImportance <- function(class_info, prune = TRUE) {
  genePairImportance = class_info$cnProc$classifier$importance
  xpairs = class_info$xpairs_list
  gene_importanceList = list()
  ignoreList = vector()
  
  # loop through individual grousp
  for(ct in names(xpairs)) {
    temp_pairList = xpairs[[ct]]
    geneNames = unique(unlist(strsplit(x = temp_pairList, split = "_"))) # get unique genes from gene pairs
    group_geneImportance = rep(0, length(geneNames))
    names(group_geneImportance) = geneNames
    
    # loop through individual gene pair
    for(tempPair in temp_pairList){
      gene1 = strsplit(tempPair, split = "_")[[1]][1]
      gene2 = strsplit(tempPair, split = "_")[[1]][2]
      
      tempImportance = genePairImportance[tempPair, 1]
      if(group_geneImportance[[gene1]] < tempImportance){
        group_geneImportance[[gene1]] = tempImportance
      }
      
      if(group_geneImportance[[gene2]] < tempImportance){
        group_geneImportance[[gene2]] = tempImportance
      }
      
    }
    
    gene_importanceList[[ct]] = group_geneImportance
    ignoreList = c(ignoreList, names(group_geneImportance))
  }
  
  # find the duplicated genes that are across different
  ignoreList = ignoreList[duplicated(ignoreList)]
  if(prune == TRUE) {
    for(ct in names(gene_importanceList)) {
      tempGene = gene_importanceList[[ct]]
      gene_importanceList[[ct]] = tempGene[!(tempGene %in% ignoreList)]
    }
  }
  
  return(gene_importanceList)
}