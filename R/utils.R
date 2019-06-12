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
extractSCE <- function(sce_object, exp_type = "normcounts"){
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
