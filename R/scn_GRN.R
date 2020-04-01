##steps for calculating GRN score and NIS score
#first reconstruct network (scn_makeGRN)
#compute for  C/T-specific GRN status (cn_barplot_grnSing)
#compute for Network influence score (NIS) (cn_nis_all)

############################################################Reconstructing network

#' @title make GRN
#' @description Construct GRN
#'
#' @param expTrain normalized expression matrix
#' @param stTrain sample table
#' @param tfs specify a list of TFs
#' @param dLevel name of the column with different categories
#' @param zThresh threshold of z score for CLR network reconstruction
#' @param computeVar boolean to indicate whether to use highly variable genes or not
#' @param dLevelGK name of the column with different germlayer categories
#' @param prune boolean limit to genes exclusively detected as CT in one CT
#' @param holm threshold of holm adjusted p value for selecting subnetwork genes
#' @param cval threshold of cval for selectin subnetwork genes. Higher cval indicates selecting higher enriched genes
#' @param method metric to construct GRN, can be pearson or MI
#' @param species when tfs is not specified, one use find tfs to specified Mm or Hs
#'
#' @return constructed GRN
#' @export
scn_makeGRN <- function(expTrain, stTrain, tfs, dLevel, zThresh = 3, 
                        computeVar=FALSE, dLevelGK = NULL, prune = FALSE, holm = 1e-4, cval=0.4, method="pearson",species="Mm"){

  # indicate target genes
  expTrain = expTrain[which(rowSums(expTrain)!=0),]
  #targetGenes = rownames(expTrain)
  
  # indicate Tfs
  if(length(tfs)==0){
    tfs=find_tfs(species)
  }
  tfs = intersect(tfs, rownames(expTrain))
  
  if(computeVar){
    cat("Computing variable genes \n")
    geneStats<-sc_statTab(expTrain, dThresh=0)
    varGenes <- findVarGenes(expTrain, geneStats = geneStats,zThresh = 2)
    cat(length(varGenes), "varGenes will be used for GRN reconstruction \n")
    expTrain = expTrain[c(varGenes,tfs),]
  }
  
  #coor_grn = grn_corr_round(expTrain)
  #zscores = grn_zscores(coor_grn, tfs)
  
  #grnall = cn_getRawGRN(zscores, coor_grn, targetGenes, zThresh=zThresh)
  grnall = scn_reconstructGRN(expDat = expTrain, tfs=tfs, method = method, zThresh = zThresh)
  
  # find preferentially expressed genes
  cat("Finished constructing general GRN", "\n")
  specGenes = cn_specGenesAll(expTrain, stTrain, holm=holm, cval=cval, dLevel=dLevel, dLevelGK = dLevelGK, prune=prune)
  
  cat("Finished finding type specific genes", "\n")
  ctGRNs = cn_specGRNs(grnall, specGenes)
  
  cat("Finished constructing type specific GRNs", "\n")
  grn_all = list(overallGRN=grnall, specGenes=specGenes,ctGRNs=ctGRNs, grnSamples=rownames(stTrain));
  
  return(grn_all)
}

#' @title reconstruct GRN
#' @description reconstruct GRN
#'
#' @param expDat normalized expression matrix
#' @param tfs specify a list of TFs
#' @param zThresh threshold of z score for CLR network reconstruction
#' @param method metric to construct GRN, can be pearson or MI
#'
#' @return reconstruct GRN
#' @export
scn_reconstructGRN <- function (expDat, tfs, method = "pearson", zThresh = 2){
  cat("Making the GRN table...\n")
  if (method == "MI") {
    ttDat = t(expDat)
    mim <- build.mim(as.data.frame(ttDat), estimator = "mi.empirical", 
                     disc = "equalwidth")
    xnet <- clr(mim)
    xcorr = cor(ttDat)
    tfsI = intersect(tfs, colnames(ttDat))
    xnet = xnet[, tfsI]
    xcorr = xcorr[, tfsI]
    grn = extractRegsDF(xnet, xcorr, rownames(expDat), zThresh)
  }
  else {
    ttDat = t(expDat)
    mim <- build.mim(ttDat, estimator = "pearson")
    xnet <- clr(mim)
    xcorr = cor(ttDat)
    tfsI = intersect(tfs, colnames(ttDat))
    xnet = xnet[, tfsI]
    xcorr = xcorr[, tfsI]
    grn =extractRegsDF(xnet, xcorr, rownames(expDat), zThresh)
  }
  
  cat("Making the iGraph...\n")
  colnames(grn)[1:2]<-c("TG", "TF")
  igTmp <-ig_tabToIgraph(grn, directed = FALSE, weights = TRUE)
  
  list(grnTable = grn, graph = igTmp)
  
}

#' find transcript factors
#'
#' find transcript factors
#' @param species defaul is 'Hs', can also be 'Mm;
#'
#' @return vector fo TF names
#'
#' @import GO.db org.Hs.eg.db
#' @export
find_tfs<-function(species='Hs'){
  
  cat("Loading gene annotations ...\n")
  require(GO.db);
  
  if(species=='Hs'){
    require(org.Hs.eg.db);
    egSymbols<-as.list(org.Hs.egSYMBOL);
    goegs<-as.list(org.Hs.egGO2ALLEGS);
  }
  
  goterms<-as.list(GOTERM);
  goids<-names(goegs);
  onts<-lapply(goids, Ontology);
  bps<-onts[onts=='BP'];
  goids<-names(unlist(bps));
  
  cat("matching gene symbols and annotations")
  gobpList<-list();
  for(goid in goids){
    egs <- goegs[[ goid ]];
    goterm<-Term(goterms[[goid]]);
    genes<-sort(unique(as.vector(unlist( egSymbols[egs] ))));
    gobpList[[goterm]]<-genes;
  }
  
  regNames<-names(gobpList)[grep("regulation of transcription", names(gobpList))];
  trs<- unique(unlist(gobpList[regNames]));
  cat("Regulation of transcription: ", length(trs),"\n");
  
  mfs<-onts[onts=='MF'];
  goidsMF<-names(unlist(mfs));
  
  gomfList<-list();
  for(goid in goidsMF){
    egs <- goegs[[ goid ]];
    goterm<-Term(goterms[[goid]]);
    genes<-sort(unique(as.vector(unlist( egSymbols[egs] ))));
    gomfList[[goterm]]<-genes;
  }
  dbs<-gomfList[['DNA binding']];
  cat("DNA binding: ", length(dbs),"\n");
  sort(intersect(trs, dbs));
}


#' @title extracts the TRs, zscores, and corr values passing thresh
#'
#' @description extracts the TRs, zscores, and corr values passing thresh
#' @param zscores zscore matrix, non-TFs already removed from columns
#' @param corrMatrix correlation matrix
#' @param genes vector of target genes
#' @param threshold zscore threshold
#'
#' @return data.frame(target=targets, reg=regulators, zscore=zscoresX, corr=correlations);
extractRegsDF<-function(zscores, corrMatrix, genes, threshold){
  
  targets<-vector();
  regulators=vector();
  zscoresX<-vector();
  correlations<-vector();
  
  targets<-rep('', 1e6);
  regulators<-rep('', 1e6);
  zscoresX<-rep(0, 1e6);
  correlations<-rep(0, 1e6);
  
  str<-1;
  stp<-1;
  for(target in genes){
    x<-zscores[target,];
    regs<-names(which(x>threshold));
    if(length(regs)>0){
      zzs<-x[regs];
      corrs<-corrMatrix[target,regs];
      ncount<-length(regs);
      stp<-str+ncount-1;
      targets[str:stp]<-rep(target, ncount);
      regulators[str:stp]<-regs;
      zscoresX[str:stp]<-zzs;
      correlations[str:stp]<-corrs;
      str<-stp+1;
    }
  }
  targets<-targets[1:stp];
  regulators<-regulators[1:stp];
  zscoresX<-zscoresX[1:stp];
  correlations<-correlations[1:stp];
  
  
  data.frame(target=targets, reg=regulators, zscore=zscoresX, corr=correlations);
}


#' @title convert a table to an igraph
#'
#' @description convert a table to an igraph. This adds an nEnts vertex attribute to count the number of entities in the sub-net
#'
#' @param grnTab table of TF, TF, maybe zscores, maybe correlations
#' @param simplify false
#' @param directed FALSE,
#' @param weights TRUE
#'
#' @return iGraph object
#' @import igraph
ig_tabToIgraph<-function(grnTab, simplify=FALSE, directed=FALSE, weights=TRUE){
  
  tmpAns<-as.matrix(grnTab[,c("TF", "TG")]);
  regs<-as.vector(unique(grnTab[,"TF"]));
  targs<-setdiff( as.vector(grnTab[,"TG"]), regs);
  
  myRegs<-rep("Regulator", length=length(regs));
  myTargs<-rep("Target", length=length(targs));
  
  types<-c(myRegs, myTargs);
  verticies<-data.frame(name=c(regs,targs), label=c(regs,targs), type=types);
  
  iG<-igraph::graph_from_data_frame(tmpAns,directed=directed,v=verticies);
  
  if(weights){
    E(iG)$weight<-grnTab$zscore;
  }
  
  if(simplify){
    iG<-simplify(iG);
  }
  V(iG)$nEnts<-1;
  iG;
}

#' @title find genes that are preferentially expressed in specified samples
#'
#' @description find genes that are preferentially expressed in specified samples
#' @param expDat expression matrix
#' @param sampTab sample table
#' @param holm sig threshold
#' @param cval R thresh
#' @param dLevel annotation level to group on
#' @param prune boolean limit to genes exclusively detected as CT in one CT
#'
#' @return a list of something
#'
#' @export
scn_findSpecGenes<-function (expDat, sampTab, holm=1e-4, cval=0.4, dLevel="description1", prune=FALSE){
  
  myPatternG<-scn_sampR_to_pattern(as.vector(sampTab[,dLevel]));
  specificSets<-apply(myPatternG, 1, cn_testPattern, expDat=expDat);
  
  # adaptively extract the best genes per lineage
  cvalT<-vector();
  ctGenes<-list();
  ctNames<-unique(as.vector(sampTab[,dLevel]));
  for(ctName in ctNames){
    x<-specificSets[[ctName]];
    
    # modification so that both upregulated and downregulated genes are going to be selected
    tmp<-rownames(x[which(x$cval>cval),]);
    tmp2<-rownames(x[which(x$holm<holm),]);
    tmp<-intersect(tmp, tmp2)
    
    up_regGenes = rep(1, length(tmp))
    names(up_regGenes) = tmp
    
    # now the downregulated genes
    tmp<-rownames(x[which(x$cval < -cval),]);
    tmp2<-rownames(x[which(x$holm<holm),]);
    tmp<-intersect(tmp, tmp2)
    
    down_regGenes = rep(-1, length(tmp))
    names(down_regGenes) = tmp
    
    totalGenes = c(up_regGenes, down_regGenes)
    
    ctGenes[[ctName]]<-totalGenes;

  }
  
  if(prune){

    specGenes<-list();
    for(ctName in ctNames){
      others<-setdiff(ctNames, ctName);
      
      exclusiveGenes<-setdiff( names(ctGenes[[ctName]]), unlist(names(ctGenes[others])));
      
      specGenes[[ctName]]<-ctGenes[[ctName]][exclusiveGenes]; # only select the exclusive genes
    }
    ans<-specGenes;
  }
  
  else{
    ans<-ctGenes;
  }
  ans;
}

#' return a pattern for use in cn_testPattern (template matching)
#'
#' return a pattern for use in cn_testPattern (template matching)
#' @param sampR vector
#'
#' @return ans
scn_sampR_to_pattern<-function (sampR){
  d_ids<-unique(as.vector(sampR));
  nnnc<-length(sampR);
  ans<-matrix(nrow=length(d_ids), ncol=nnnc);
  for(i in seq(length(d_ids))){
    x<-rep(0,nnnc);
    x[which(sampR==d_ids[i])]<-1;
    ans[i,]<-x;
  }
  colnames(ans)<-as.vector(sampR);
  rownames(ans)<-d_ids;
  ans;
}

#' @title template matching
#'
#' @description test correlation between idealized expression pattern and target gene
#' @param pattern vector of pattern
#' @param expDat expression matrix
#'
#' @return data.frame(row.names=geneids, pval=pval, cval=cval,holm=holm);
cn_testPattern<-function(pattern, expDat){
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

#' finds general and context dependent specifc genes
#'
#' finds general and context dependent specifc genes
#' @param expDat expression matrix
#' @param sampTab sample table
#' @param holm pvalue threshold for template matching
#' @param cval template matching threshold for overall CT specific expression
#' @param cvalGK template matching threshold for developmentally shared CT specific expression
#' @param dLevel "description1",
#' @param dLevelGK "description2"
#'
#' @return list of $matcher${cell_type}->{germ_layer}$context$general${cell_type}->gene vector etc
#' @export
cn_specGenesAll<-function(expDat, sampTab,holm=1e-4,cval=0.4,cvalGK=0.75, dLevel, dLevelGK=NULL,prune=FALSE){
  matcher<-list();
  general<-scn_findSpecGenes(expDat, sampTab, holm=holm, cval=cval, dLevel=dLevel,prune=prune);
  ctXs<-list()# one per germlayer
  if(!is.null(dLevelGK)){
    
    germLayers<-unique(as.vector(sampTab[,dLevelGK]));
    for(germlayer in germLayers){
      stTmp<-sampTab[sampTab[,dLevelGK]==germlayer,];
      expTmp<-expDat[,rownames(stTmp)];
      xxx<-scn_findSpecGenes(expTmp, stTmp, holm=holm, cval=cvalGK,dLevel=dLevel, prune=prune);
      cts<-names(xxx);
      for(ct in cts){
        matcher[[ct]]<-germlayer;
        # remove general ct-specific genes from this set
        a<-general[[ct]];
        b<-xxx[[ct]];
        ba<-setdiff(b, a);
        both<-union(a,b);
        xxx[[ct]]<-ba;
      }
      ctXs[[germlayer]]<-xxx;
    }
  }
  ctXs[['general']]<-general;
  list(context=ctXs, matcher=matcher);
}

#' extract sub-networks made up of CT genes;
#'
#' extract sub-networks made up of CT genes;
#' @param rawGRNs result of running cn_getRawGRN
#' @param specGenes result of running cn_specGenesAll
#'
#' @return list(geneLists=geneLists, graphLists=graphLists, tfTargets=tfTargets)
#' @export
cn_specGRNs<-function(rawGRNs, specGenes){
  
  # should return a list of gene lists and igraphs
  geneLists<-list();
  graphLists<-list();
  
  groupNames<-names(specGenes[['context']][['general']]);
  
  big_graph<-rawGRNs[['graph']];
  
  matcher<-specGenes$matcher;
  
  allgenes<-V(big_graph)$name;
  
  for(ct in groupNames){
    cat(ct,"\n")
    if(!is.null(names(matcher))){ #TODO modifiy the germ line stuff later
      gll<-matcher[[ct]];
      cat(ct," ",gll,"\n");
      mygenes<-union(specGenes[['context']][['general']][[ct]], specGenes[['context']][[gll]][[ct]]);
    }
    else{
      mygenes<-names(specGenes[['context']][['general']][[ct]])
    }
    
    geneLists[[ct]]<-specGenes[['context']][['general']][[ct]][intersect(allgenes, mygenes)];
    
    graphLists[[ct]]<-induced.subgraph(big_graph, names(geneLists[[ct]]));
    
  }
  
  tfTargets<-makeTLs(graphLists);
  
  list(geneLists=geneLists, graphLists=graphLists, tfTargets=tfTargets);
}

#' get targets of tFs
#'
#' get targets of tFs
#' @param graphList a list of networks represented as iGraphs
#'
#' @return list of tf=>targets
makeTLs<-function(graphList){
  
  tfTargs<-list();
  nnames<-names(graphList);
  for(nname in nnames){
    tfTargs[[nname]]<-get_targets(graphList[[nname]]); # get the name of the genelist
  }
  tfTargs;
}

#' get targets of a tf
#'
#' get targets of a tf
#' @param aGraph an iGraph
#'
#' @return target list
get_targets<-function(aGraph){
  
  targList<-list();
  regs<-V(aGraph)$label[V(aGraph)$type=='Regulator'];
  if(length(regs)>0){
    for(reg in regs){
      targList[[reg]]<-unique(sort(V(aGraph)$label[neighbors(aGraph, reg)]));
    }
  }
  targList;
}
