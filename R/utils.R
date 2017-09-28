# scCellNet
# (C) Patrick Cahan 2012-2017

# commonly used or misc functions

#' @export
getGenesFromGO<-function# return the entrez gene ids of a given a GOID, for now assumes mouse
(GOID, # GO id to find genes for
 annList 
){
  sort(as.vector(unlist(annList[['egSymbols']][annList[['goegs']][[GOID]]])));
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


#' find transcript factors
#'
#' find transcript factors
#' @param annotation
#' @param species defaul is 'Hs', can also be 'Mm;
#' @param ontology default is BP 
#'
#' @return vector fo TF names
#' @export
#' @importFrom AnnotationDbi as.list
#'
find_genes_byGo<-function# 
(annotation,
  species='Hs',
  onto="BP"
){

  cat("Loading gene annotations ...\n")
  require(GO.db);

  if(species=='Hs'){
    require(org.Hs.eg.db);
    egSymbols<-as.list(org.Hs.egSYMBOL);
    goegs<-as.list(org.Hs.egGO2ALLEGS);
  }
  else{
    require(org.Mm.eg.db);
    egSymbols<-as.list(org.Mm.egSYMBOL);
    goegs<-as.list(org.Mm.egGO2ALLEGS);
  }

  goterms<-as.list(GOTERM);
  goids<-names(goegs);
  onts<-lapply(goids, Ontology);
  bps<-onts[onts==onto];
  goids<-names(unlist(bps));

  cat("matching gene symbols and annotations")
  gobpList<-list();
  for(goid in goids){
    egs <- goegs[[ goid ]];
    goterm<-Term(goterms[[goid]]);
    genes<-sort(unique(as.vector(unlist( egSymbols[egs] ))));
    gobpList[[goterm]]<-genes;
  }

  ### newHsTRs<-gobpList[['regulation of transcription, DNA-dependent']];
  regNames<-names(gobpList)[grep(annotation, names(gobpList))];
  trs<- unique(unlist(gobpList[regNames]));
  cat(annotation, ": ", length(trs),"\n");
  sort(trs)

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
