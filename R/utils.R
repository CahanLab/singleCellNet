# scCellNet
# (C) Patrick Cahan 2012-2017

# commonly used or misc functions

#' @export
ctMerge<-function(sampTab, annCol, ctVect, newName){
  oldann<-as.vector(sampTab[,annCol])
  newann<-oldann
  for(oldName in ctVect){
    xi<-which(oldann==oldName)
    newann[xi]<-newName
  }
  xnot<-which(colnames(sampTab)!=annCol)
  ans<-sampTab[,xnot]
  cbind(ans, newAnn=newann)
}

#' @export
ctRename<-function(sampTab, annCol, oldName, newName){
  oldann<-as.vector(sampTab[,annCol])
  newann<-oldann
  xi<-which(oldann==oldName)
  newann[xi]<-newName
  xnot<-which(colnames(sampTab)!=annCol)
  ans<-sampTab[,xnot]
  cbind(ans, newAnn=newann)
} 

#' @export
splitCommon<-function(sampTab, ncells, dLevel="cell_ontology_class"){
  cts<-unique(as.vector(sampTab[,dLevel]))
  trainingids<-vector()
  for(ct in cts){
    cat(ct,": ")
    stX<-sampTab[sampTab[,dLevel]==ct,]
    ccount<-nrow(stX)-3
    ccount<-min(ccount, ncells)
    cat(nrow(stX),"\n")
    trainingids<-append(trainingids, sample(rownames(stX), ccount))
  }
  val_ids<-setdiff(rownames(sampTab), trainingids)
  list(train=sampTab[trainingids,], val=sampTab[val_ids,])
}


#' @export
loadLoomExp<-function# load a loom object containing expression data
(path,
  cellNameCol='obs_names'
  ){

  lfile <- connect(filename = path)
  geneNames<-lfile[["row_attrs"]][["var_names"]][]
  cellNames<-lfile[["col_attrs"]][["obs_names"]][]
  expMat<- t(lfile[["matrix"]][1:length(cellNames),])
  rownames(expMat)<-geneNames
  colnames(expMat)<-cellNames
  expMat
}

#' @export
loadLoomExpCluster<-function# load a loom object containing expression  + cluster info
(path,
  cellNameCol='obs_names',
  xname='cluster'
  ){
  lfile <- connect(filename = path, skip.validate = TRUE)
  geneNames<-lfile[["row_attrs"]][["var_names"]][]
  cellNames<-lfile[["col_attrs"]][["obs_names"]][]
  expMat<- t(lfile[["matrix"]][1:length(cellNames),])
  rownames(expMat)<-geneNames
  colnames(expMat)<-cellNames
  
  cluster_old <- lfile[['col_attrs']][[xname]][]

  sampTab <- data.frame(cell_name=cellNames, cluster=cluster_old)
  row.names(sampTab) <- cellNames
  lfile$close_all()
  list(expDat = expMat, sampTab = sampTab)
}




#' @export
getGenesFromGO<-function# return the entrez gene ids of a given a GOID, for now assumes mouse
(GOID, # GO id to find genes for
 annList 
){
  sort(as.vector(unlist(annList[['egSymbols']][annList[['goegs']][[GOID]]])));
}



dumbfunc<-function(
 aNamedList){
  ans<-vector()
  nnames<-names(aNamedList)
  for(nname in nnames){
    ans<-append(ans, rep(nname, length(aNamedList[[nname]])))
  }
  ans
}

#' row average (or median) based on groups
#'
#' row average (or median) based on groups
#' @param exp expression df
#' @param groupings groupings
#' @param type mean or media
#'
#' @return return a dataframe of mean or median-ed data based on given groupings.  colnames become the column name of the first sample in each group from the original data
#'
#' @export
GEP_makeMean<-function
(exp,
 groupings,
 type='mean'
){
  
  
  ans<-data.frame();
  grps<-unique(groupings);
  if(type=='mean'){
    for(grp in grps){
      gi<-which(groupings==grp);
      if(length(gi)==1){
        
        if(nrow(ans)==0){
          ans<-data.frame(exp[,gi]);
        }else{
          ans<-cbind(ans, exp[,gi]);
        }
      }
      else{
        xxx<-apply(exp[,gi],1,mean);
        if(nrow(ans)==0){
          ans<-data.frame(xxx);
        }
        else{
          ans<-cbind(ans, xxx);
        }
      }
    }
  }
  else{
    for(grp in grps){
      gi<-which(groupings==grp);
      xxx<-apply(exp[,gi],1,median);
      if(nrow(ans)==0){
        ans<-data.frame(xxx);
      }
      else{
        ans<-cbind(ans, xxx);
      }
    }
  }
  
  colnames(ans)<-grps;
  ans;
  ### data.frame of mean or median-ed data based on given groupings
}



# get GO:IDs
#' @export
annSetUp<-function(){
library("org.Mm.eg.db")
library(GO.db)
  cat("Getting GO:IDs\n");
  goegs<-as.list(org.Mm.egGO2ALLEGS);

  # get EG <-> Symbols
  cat("Getting Symbols\n");
  egSymbols<-as.list(org.Mm.egSYMBOL);
  list(goegs=goegs, egSymbols=egSymbols)
}



#' 1-PCC distance
#'
#' 1-PCC distance
#' @param x numeric matrix
#' 
#' @return distance matrix  
#'
#' @examples
#' xdist<-utils_myDist(t(expDat))
#' plot(hclust(xdist, 'ave'), hang=-1)
#'
#' @export
utils_myDist<-function
(x
){
  as.dist(1-cor(t(x)));
}

#' loads an R object when you don't know the name
#'
#' loads an R object when you don't know the name
#' @param fname file
#'
#' @return variable
#'
#' @export
utils_loadObject<-function
(fname
 ### file name
){
  x<-load(fname);
  get(x);
}

#' strip whitespace from a string
#'
#' strip whitespace from a string
#' @param string string
#'
#' @return new string
#'
#' @export
utils_stripwhite<-function
### 
(string
 #### string
 ){
  gsub("^\\s+|\\s+$", "", string)
}

#' print date
#'
#' print date
#' @return string
#'
#' @export
utils_myDate<-function
### 
()
{
  format(Sys.time(), "%b_%d_%Y");
}

#' reduces full path to filename
#'
#' reduces full path to filename
#' @param string
#'
#' @return something
#'
#' @export
utils_strip_fname<-function #
(str){
  a<-strsplit(str, "/")[[1]];
  a[length(a)];
}

utils_stderr<-function
### calculate standard error
(x){
  sqrt(var(x)/length(x));
  ### stderr
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


zscoreVect<-function
### Compute the mean zscore of given genes in each sample
(genes,
 ### genes
 xvals,
 ### named vector
 tVals,
 ### tvals
 ctt
 ### ctt
 ){
  ans<-vector();
  for(gene in genes){
    ans<-append(ans, zscore(xvals[gene], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]]));
  }
  ans;
  ### zscore vector
}

#' make Inf and -Inf values sensible
#'
#' make Inf and -Inf values sensible
#' @param zMat zMat
#'
#' @return corrected zMat
#'
#' @export
cn_correctZmat<-function
(zmat){
  myfuncInf<-function(vect){
    xi<-which(vect=='Inf')
    if(any(xi)){
      mymax<-max(vect[-xi])
      vect[xi]<-mymax
    }
    vect
  }
  zmat<-apply(zmat,2, myfuncInf)
  zmat[is.na(zmat)]<-0
  zmat
}

#' weighted subtraction from mapped reades
#'
#' Simulate expression profile of  _total_ mapped reads
#' @param vector of total mapped reads per gene/transcript
#' @param total post transformation sum of read counts
#'
#' @return vector of downsampled read mapped to genes/transcripts
#'
#' @export
downSampleW<-function
(vector,
 total=1e5,
 dThresh=0){ 

  totalSignal<-sum(vector)
  wAve<-vector/totalSignal
###  resid<-sum(vector)-total #num to subtract from sample
  resid<-totalSignal-total #num to subtract from sample
  residW<-wAve*resid # amount to substract from each gene
  ans<-vector-residW
  ans[which(ans<dThresh)]<-0
  ans
}

#' weighted subtraction from mapped reades, applied to all
#'
#' Simulate expression profile of  _total_ mapped reads
#' @param expRaw matrix of total mapped reads per gene/transcript
#' @param total numeric post transformation sum of read counts
#'
#' @return vector of downsampled read mapped to genes/transcripts
#'
#' @export
weighted_down<-function
(expDat,
 total,
 dThresh=0
 ){
  if(class(expDat)[1]!='matrix'){
    cSums  <- Matrix::colSums(expDat)
    props <- Matrix::t(expDat) / cSums
    rrids  <- cSums - total
    tmpAns <- expDat - Matrix::t(props * rrids)
    tmpAns[Matrix::which(tmpAns<dThresh)] <- 0
  }
  else{
    cSums  <- colSums(expDat)
    props <- t(expDat) / cSums
    rrids  <- cSums - total
    tmpAns <- expDat - t(props * rrids)
    tmpAns[which(tmpAns<dThresh)] <- 0
  }
  
  tmpAns
}


#' weighted subtraction from mapped reades and log transform the data, applied to all
#'
#' Simulate expression profile of  _total_ mapped reads
#' @param expRaw matrix of total mapped reads per gene/transcript
#' @param total numeric post transformation sum of read counts
#'
#' @return vector of downsampled read mapped to genes/transcripts
#'
#' @export
trans_prop<-function
(expDat,
 total,
 dThresh=0
 ){
  if(class(expDat)[1]!='matrix'){
    cSums  <- Matrix::colSums(expDat)
    props <- Matrix::t(expDat) / cSums
    rrids  <- cSums - total
    tmpAns <- expDat - Matrix::t(props * rrids)
    tmpAns[Matrix::which(tmpAns<dThresh)] <- 0
  }
  else{
    cSums  <- colSums(expDat)
    props <- t(expDat) / cSums
    rrids  <- cSums - total
    tmpAns <- expDat - t(props * rrids)
    tmpAns[which(tmpAns<dThresh)] <- 0
  }
  ans = log(1+tmpAns)
  ans
}

#' @export
trans_zscore<-function
(expRaw){
  apply(expRaw, 2, scale)
}

#' @export
trans_binarize<-function
(expRaw,
  threshold=1){
  expRaw[expRaw<threshold]<-0
  expRaw[expRaw>0]<-1
  expRaw
}

#' extract sampTab and expDat seurat object into regular S3 objects
#' @param seurat_object
#' @param exp_slot_name
#' @return list
#' @export
extractSeurat <- function(seurat_object, exp_slot_name = "counts"){
  
  #extract metadata
  sampTab = seurat_object@meta.data
  
  #extract expression matrix
  expDat = as.matrix(GetAssayData(seurat_object, slot = exp_slot_name))   
  
  return(list(sampTab = sampTab, expDat = expDat))
  
}

#' extract sampTab and expDat sce object into regular S3 objects
#' @param sce_object
#' @param exp_type
#' @param list
#' @export
extractSCE <- function(sce_object, exp_slot_name = "counts"){
  #extract metadata
  sampTab = as.data.frame(colData(sce_object, internal = TRUE))
  sampTab$sample_name = rownames(sampTab)
  
  #extract expression matrix
  if(exp_type == "counts"){
    expDat = counts(sce_object)
  }
  
  if(exp_type == "normcounts"){
    expDat = normcounts(sce_object)
  }
  
  if(exp_type == "logcounts"){
    expDat = logcounts(sce_object)
  }
  
  return(list(sampTab = sampTab, expDat = expDat))
  
}

#' @export
getUniqueGenes<-function(genes, transID='id', geneID='symbol'){ #
  rownames(genes)<-as.vector(genes[,transID])
  dgenes<-names(which(table(genes[,2])>1))
  xipos<-which(genes[,geneID] %in% dgenes)
  ###xipos<-match(dgenes, genes[,geneID])
  tnames<-rownames(genes)
  sgenes<-setdiff(tnames, tnames[xipos])
  genes2<-genes[sgenes,]
  xx<-tnames[match(dgenes, genes[,geneID])]
  genes2<-rbind(genes2, genes[xx,])
  genes2
}

#' @export
removeRed<-function(expOb,transID="id", geneID="symbol"){
     genes<-getUniqueGenes( fData(expOb), transID, geneID)
     expOb[rownames(genes),]
}

#' @export
load10x_mtx<-function
(path,
 mtx_fname="matrix.mtx",
 gene_fname="genes.tsv",
 bc_fname="barcodes.tsv",
 bcFiltered=TRUE
 ){
        mat_fn<-paste0(path,"/",mtx_fname)
        gene_fn<-paste0(path,"/",gene_fname)
        barcode_fn<-paste0(path,"/",bc_fname)
        summary_fn<-NULL
        res <- load_cellranger_matrix_from_files(mat_fn, gene_fn, barcode_fn, summary_fn)
        res@barcode_filtered <- bcFiltered
    res@subsampled <- FALSE
    #res@pipestance_path <- pipestance_path
    res
}

#' @export
load10x_mtx<-function
(path,
 mtx_fname="matrix.mtx",
 gene_fname="genes.tsv",
 bc_fname="barcodes.tsv",
 bcFiltered=TRUE
 ){
 	mat_fn<-paste0(path,"/",mtx_fname)
 	gene_fn<-paste0(path,"/",gene_fname)
 	barcode_fn<-paste0(path,"/",bc_fname)
 	summary_fn<-NULL
	res <- load_cellranger_matrix_from_files(mat_fn, gene_fn, barcode_fn, summary_fn)
	res@barcode_filtered <- bcFiltered
    res@subsampled <- FALSE
    #res@pipestance_path <- pipestance_path
    res
}

#' Loads 10x data mtx
#'
#' @param path path
#' @param prefix label to facilitate merging data sets when UMIs might collide
#' @param mtx_fname mtx_fname
#' @param gene_fname gene_fname
#' @param bc_fname bc_fname
#' @param bcFiltered bcFiltered
#' @param removeRedundant removeRedundant
#'
#' @return list of list(sampTab=sampTab, expDat=expDat)
#'
#' @export
load10x<-function
(path,
 prefix,
 mtx_fname="matrix.mtx",
 gene_fname="genes.tsv",
 bc_fname="barcodes.tsv",
 bcFiltered=TRUE,
 removeRedundant=TRUE
 ){
	bigDat<-load10x_mtx(path, mtx_fname, gene_fname, bc_fname, bcFiltered)
	redDat<-removeRed(bigDat)

	geneTab<-fData(redDat)
	rownames(geneTab)<-geneTab[,2]

	expDat<-as.matrix(exprs(redDat))
	rownames(expDat)<-rownames(geneTab)

	# sample_name == prefix + inc
	# sample_id == prefix + inc + barcode
	
	# barcode == barcode
	# colnames of expDat == sample_name
	# rownames(sampTab) = sample_name

	snames<-paste(prefix, 1:ncol(expDat), sep="_")
	sids<-paste(snames, colnames(expDat), sep="_")

	umis<-apply(expDat, 2, sum)
	sampTab<-data.frame(sample_id=sids, sample_name=snames, barcode=colnames(expDat), umis=umis, prefix=rep(prefix, length(umis)))
	rownames(sampTab)<-snames
	colnames(expDat)<-snames

	list(sampTab=sampTab, expDat=expDat)
}

#' Loads 10x data mtx
#'
#' @param basePath path
#' @param pNames are prefix labels to facilitate merging data sets when UMIs might collide
#' @param nCells number of cells per mtx to load
#' @param secPath secPath
#'
#' @return list of list(sampTab=sampTab, expDat=expDat)
#'
#' @export
mergeLoad10x<-function
(basePath,
 pNames,
 nCells=1e3,
 secPath="filtered_matrices_mex/hg19"){

	expList<-list()
	for(pName in pNames){
		cat(pName,"\n")
		ppath<-paste0(basePath,pName,"/",secPath)
		tmpX<-load10x(ppath, pName)
		if(nCells>0){
			cells<-sample(rownames(tmpX[['sampTab']]), nCells)
		}
		else{
			cells<-rownames(tmpX[['sampTab']])
		}
		expList[[pName]]<- list(sampTab=tmpX[['sampTab']][cells,], expDat=tmpX[['expDat']][,cells])
		cat(nrow(expList[[pName]][['expDat']]), "\n")
	}
	cgenes<-rownames(expList[[1]][['expDat']])
	for(pName in pNames){
		cgenes<-intersect(cgenes, rownames(expList[[pName]][['expDat']]))
	}
	cat("Number of genes:",length(cgenes),"\n")

	if(nCells>0){
		nCellsTotal<-nCells*length(pNames)
	}
	else{ # load them all
		nCellsTotal<- sum( unlist(lapply(expList, function(alist){ ncol(alist[['expDat']])})))
	}
	cat("Number of cells: ",nCellsTotal,"\n")

	expAll<-matrix(0, nrow=length(cgenes), ncol=nCellsTotal)
	rownames(expAll)<-cgenes
	stAll<-data.frame()
	str<-1
	##stp<-nCells
	cat(pNames[1],"\n")
	cat("named:",ncol(expList[[pNames[1]]][['expDat']]),"\n")
	cat("numbered:",ncol(expList[[1]]),"\n")	
	stp<-ncol(expList[[pNames[1]]][['expDat']])
#	for(pName in pNames){
	for(i in seq(length(pNames))){
		pName<-pNames[i]
		cat(str,"-",stp,"\n")
		expAll[cgenes,str:stp]<-expList[[pName]][['expDat']][cgenes,]
		stAll<-rbind(stAll, expList[[pName]][['sampTab']])
		str<-stp+1

		##stp<-str+nCells-1
		if(i<length(pNames)){
			inc<-ncol(expList[[ pNames[i+1] ]][['expDat']])
		}
		stp<-str+inc-1
	}
	colnames(expAll)<-rownames(stAll)
	list(sampTab=stAll, expDat=expAll)

}


