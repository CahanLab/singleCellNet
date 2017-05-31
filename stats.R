# Patrick Cahan (C) 2017
# patrick.cahan@gmail.com


findDepth<-function# count 'reads across a matrix
(aMat){
  apply(aMat, 2, sum)
}

rankTrans<-function# rank genes,  log
(expDat){
  rranks<-apply(expDat,2,rank, ties='last')
  #rranks<-maxR-rranks
  log(1+rranks)
}


ocells<-function # order cells within a group by a variable
(sampTab,
 expX,
 oList, # grpname=>list(gene, oType)
 cname="group"
){

 ans<-vector()
 for(grp in names(oList)){
    sids<-rownames(sampTab[which(sampTab[,cname]==grp),])
    cat(length(sids),"\n")
    otype<-oList[[grp]][['oType']]
    gene<-oList[[grp]][['gene']]
    cat(gene, " ", otype,"\n")
    tmpans<-colnames(expX[,sids])[(order(expX[gene,sids], decreasing=otype))]
    cat(length(tmpans),"\n")
    ans<-append(ans, tmpans)
  }
  ans
}




minTab<-function #subsample the table
(sampTab, dLevel){
  myMin<-min(table(sampTab[,dLevel]))
  nST<-data.frame()
  grps<-unique(as.vector(sampTab[,dLevel]))
  for(grp in grps){
    stX<-sampTab[sampTab[,dLevel]==grp,]
    nST<-rbind(nST, stX[sample(rownames(stX), myMin),])
  }
  nST
}

sc_sampR_to_pattern<-function#
(sampR){
  d_ids<-unique(as.vector(sampR));
  nnnc<-length(sampR);
#  ans<-matrix(nrow=length(d_ids), ncol=nnnc);
  ans<-list()
  for(d_id in d_ids){
    x<-rep(0,nnnc);
    x[which(sampR==d_id)]<-1;
    ans[[d_id]]<-x;
  }
  ans
}


sc_testPattern<-function(pattern, expDat){
  pval<-vector();
  cval<-vector();
  geneids<-rownames(expDat);
  llfit<-ls.print(lsfit(pattern, t(expDat)), digits=25, print=FALSE);
  xxx<-matrix( unlist(llfit$coef), ncol=8,byrow=TRUE);
  ccorr<-xxx[,6];
  cval<- sqrt(as.numeric(llfit$summary[,2])) * sign(ccorr);
  pval<-as.numeric(xxx[,8]);

  #qval<-qvalue(pval)$qval;
  holm<-p.adjust(pval, method='holm');
  #data.frame(row.names=geneids, pval=pval, cval=cval, qval=qval, holm=holm);
  data.frame(row.names=geneids, pval=pval, cval=cval,holm=holm);
}

getTopGenes<-function
(xdat,
  topN=3)
# cval=.35,
# holm=1e-2)
{
#  qqq<-xdat[which(xdat$cval>cval & xdat$holm<holm),]
#  rownames(qqq[order(qqq$cval, decreasing=TRUE),][1:topN,])
  rownames(xdat[order(xdat$cval, decreasing=TRUE),][1:topN,])
}


par_findSpecGenes<-function#
(expDat, ### expression matrix
  sampTab, ### sample tableß
  ###holm=1e-50, ### sig threshold
  ###cval=0.5, ### R thresh
  dLevel="group", #### annotation level to group on
  prune=FALSE, ### limit to genes exclusively detected as CT in one CT
  minSet=TRUE,
  ncore=4
){

  newST<-sampTab
  if(minSet){
    cat("Making reduced sampleTable\n")
    newST<-minTab(sampTab, dLevel)
  }
  
  ans<-list()

  cat("Making patterns\n")
  myPatternG<-sc_sampR_to_pattern(as.vector(newST[,dLevel]));
  expDat<-expDat[,rownames(newST)]
  
  cat("Testing patterns\n")
 # aClust<-parallel::makeCluster(ncore, type='FORK')
  specificSets<-lapply(myPatternG, sc_testPattern, expDat=expDat)
 # stopCluster(aClust)
  cat("Done testing\n")
  
  

  if(FALSE){
  # adaptively extract the best genes per lineage
  cvalT<-vector();
  ctGenes<-list();
  ctNames<-unique(as.vector(sampTab[,dLevel]));
  for(ctName in ctNames){
    x<-specificSets[[ctName]];
#    tmp<-rownames(x[x$cval>cval,]);
#    tmp2<-rownames(x[x$holm<holm,]);
     qqq<-x[which(x$cval>cval & x$holm<holm),]
     x<-qqq[order(qqq$cval, decreasing=TRUE),]
     tmp<-rownames(x)
     
#    tmp<-intersect(tmp, tmp2)
    ctGenes[[ctName]]<-tmp;
###    cvalT<-append(cvalT, cval);
  }

  if(prune){
  # now limit to genes exclusive to each list
   specGenes<-list();
   for(ctName in ctNames){
     others<-setdiff(ctNames, ctName);
     x<-setdiff( ctGenes[[ctName]], unlist(ctGenes[others]));
     specGenes[[ctName]]<-x;
   }
   ans<-specGenes;
  }
  else{
   ans<-ctGenes;
  }
  names(ans)<-as.character(ctNames)
  }
  specificSets
}




summModBin<-function#
(expMat){
  apply(expMat, 2, mean)
}

summModsBin<-function(expMat, geneLists, thresh=0){
  expBin<-expMat
  expBin[expBin>thresh]<-1

  ans<-matrix(0,nrow=length(geneLists), ncol=ncol(expMat))
  for(i in seq(length(geneLists))){
    ans[i,]<-summModBin(expBin[geneLists[[i]],])
  }
  rownames(ans)<-names(geneLists)
  colnames(ans)<-colnames(expMat)
  ans
}


summByzmean<-function # compute summary value based on mean of zscores
(expMat){
  zscores<-t(scale(t(expMat)))
  apply(zscores, 2, mean)
}


summMeanZscores<-function(expMat, geneLists){
  ans<-matrix(0,nrow=length(geneLists), ncol=ncol(expMat))
  for(i in seq(length(geneLists))){
    ans[i,]<-summByzmean(expMat[geneLists[[i]],])
  }
  rownames(ans)<-names(geneLists)
  colnames(ans)<-colnames(expMat)
  ans
}


summMod<-function#
(expMat){
  myfunc<-function(avect, maxes){ sum(avect/maxes) } 
  maxes<-apply(expMat, 1, max)
  apply(expMat,2, myfunc, maxes)
  
}

summMods<-function(expMat, geneLists){
  ans<-matrix(0,nrow=length(geneLists), ncol=ncol(expMat))
  for(i in seq(length(geneLists))){
    ans[i,]<-summMod(expMat[geneLists[[i]],])
  }
  rownames(ans)<-names(geneLists)
  colnames(ans)<-colnames(expMat)
  ans
}

sc_statTab<-function# make a gene stats table (i.e. alpha, mu, etc)
(expDat, # expression matrix
 dThresh=0 # threshold for detection
 ){
  
  statTab<-data.frame()
  muAll<-sc_compMu(expDat, threshold=dThresh);
  alphaAll<-sc_compAlpha(expDat,threshold=dThresh);
  meanAll<-apply(expDat, 1, mean);
  covAll<-apply(expDat, 1, sc_cov);
  fanoAll<-apply(expDat,1, sc_fano);
  maxAll<-apply(expDat, 1, max);
  sdAll<-apply(expDat, 1, sd);
  
  statTabAll<-data.frame(gene=rownames(expDat), mu=muAll, alpha=alphaAll, overall_mean=meanAll, cov=covAll, fano=fanoAll, max_val=maxAll, sd=sdAll)
  statTabAll;
}


sc_compAlpha<-function # compute alpha given detection threshold
(expMat,
 threshold=0,
  pseudo=FALSE){
  
  indexFunction<-function(vector, threshold){
    names(which(vector>threshold));
  }
  
  indexes<-apply(expMat, 1, indexFunction, threshold);
  alphas<-unlist(lapply(indexes, length));
  ans<-alphas/ncol(expMat)
  if(pseudo){
    ans<-(alphas+1)/(ncol(expMat)+1)
  }
  ans
}

sc_compMu<-function # compute Mu given threshold
(expMat,
 threshold=0){
  
  afunct<-function(vector, threshold){
    mean( vector[which(vector>threshold)] );
  } 
  
  mus<-unlist(apply(expMat, 1, afunct, threshold))
  mus[is.na(mus)]<-0;
  mus;
}

repNA<-function# replavce NAs with 0
(vector){
  vector[which(is.na(vector))]<-0;
  vector;
}

sc_fano<-function #compute fano factor on vector
(vector){
  var(vector)/mean(vector);
}

sc_cov<-function# compute coeef of variation on vector
(vector){
  sd(vector)/mean(vector);
}


#' find genes that pass criteria
#'
#' based on idea that reliably detected genes will either be detected in many cells, or highly expressed in a small cels of cells (or both
#'
#' @param geneStats result of running sc_statTab
#' @param alpha1 proportion of cells in which a gene must be considered detected (as defined in geneStats)
#' @param alpha2 lower proportion of cells for genes that must have higher expression level
#' @param mu threshold, average expression level of genes passing the lower proportion criteria
#'
#' @return vector of gene symbols
#' 
#' @export
#'
sc_filterGenes<-function
(geneStats,
 alpha1=0.1,
 alpha2=0.01,
 mu=2){
	passing1<-rownames(geneStats[geneStats$alpha>alpha1,])
	notPassing<-setdiff(rownames(geneStats), passing1)
	geneStats<-geneStats[notPassing,]
	c(passing1, rownames(geneStats[which(geneStats$alpha>alpha2 & geneStats$mu>mu),]))
}


#' find cells that pass criteria
#'
#' based purely on umis
#'
#' @param sampTab, which must have UMI column
#' @param minVal umis must exceed this
#' @param maxValQuant quantile to select max threshold
#'
#' @return vector rownames(sampTab) meeting criteria
#' 
#' @export
#'
sc_filterCells<-function
(sampTab,
 minVal=1e3,
 maxValQuant=0.95){
	stX<-sampTab[sampTab$umis>minVal,]
	qThresh<-quantile(sampTab$umis, maxValQuant)
	rownames(stX[stX$umis<qThresh,])
}



#' finds genes higher in one group vs others
#'
#' finds genes higher in one group vs others
#'
#' @param expDat expression matrix
#' @param list of cells in eahc groupsampTab result of running cutreeDynamicTree 
#' @param 
#'
#' @return list of grps -> names
#' 
#' @export
#'
sc_findEnr<-function #
(expDat,
 sampTab,
 dLevel="group")
{
  groups<-unique(as.vector(sampTab[,dLevel]))
  myMeans<-matrix(0,nrow=nrow(expDat), ncol=length(groups))
  rownames(myMeans)<-rownames(expDat)
  colnames(myMeans)<-groups

  for(group in groups){
    xi<-which(sampTab[,dLevel]==group)
    myMeans[,group]<-apply(expDat[,xi], 1, median)
  }

  ans<-list()
  for(group in groups){
    cat(group,"\n")
    others<-setdiff(groups, group)
    tmpMeans<-apply(myMeans[,others], 1, median)
    myDiff<-myMeans[,group] - tmpMeans
    ans[[group]]<-rownames(expDat)[order(myDiff, decreasing=TRUE)]
  }
  ans;
}

binGenesAlpha<-function### bins genes into x groups based on overallmean
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

binGenes<-function### bins genes into x groups based on overallmean
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
    xnames<-rownames(geneStats[which(geneStats$overall_mean<=borders[i]),]);
    binGroup[xnames]<-i;
  }
  cbind(geneStats, bin=binGroup);
}

# find most variable genes
findVarGenes<-function(expNorm, geneStats,zThresh=2, meanType="overall_mean"){
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
  sg<-binGenes(geneStats, meanType=meanType)
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



