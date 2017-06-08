# Patrick Cahan (C) 2017
# patrick.cahan@gmail.com


#' @export
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


# compute alpha given detection threshold
#' @export
sc_compAlpha<-function 
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

# compute Mu given threshold
#' @export
sc_compMu<-function 
(expMat,
 threshold=0){
  
  afunct<-function(vector, threshold){
    mean( vector[which(vector>threshold)] );
  } 
  
  mus<-unlist(apply(expMat, 1, afunct, threshold))
  mus[is.na(mus)]<-0;
  mus;
}

# replavce NAs with 0
repNA<-function
(vector){
  vector[which(is.na(vector))]<-0;
  vector;
}

#compute fano factor on vector
sc_fano<-function
(vector){
  var(vector)/mean(vector);
}

# compute coeef of variation on vector
sc_cov<-function
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
#' @param dLevel dLevel
#'
#' @return list of grps -> names
#' 
#' @export
#'
sc_findEnr<-function 
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

#' @export
enrDiff<-function 
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

  ans<-matrix(0,nrow=nrow(expDat), ncol=length(groups))
  colnames(ans)<-groups
  rownames(ans)<-rownames(expDat)
  for(group in groups){
    cat(group,"\n")
    others<-setdiff(groups, group)
    tmpMeans<-apply(myMeans[,others], 1, median)
    myDiff<-myMeans[,group] - tmpMeans
    ans[,group]<-myDiff
  }
  ans
}



### bins genes into x groups based on overallmean
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
    xnames<-rownames(geneStats[which(geneStats$overall_mean<=borders[i]),]);
    binGroup[xnames]<-i;
  }
  cbind(geneStats, bin=binGroup);
}

# find most variable genes
#' @export
findVarGenes<-function
(expNorm,
  geneStats,
  zThresh=2, 
  meanType="overall_mean"){
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

