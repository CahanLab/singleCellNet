#pass the clustering step to the classification
#creating a classifier at specific level
#using differential expression to build the classifiers
#only build classifiers for the top and bottom level 

pipe_butter_gpa_diffExp_cn_make <-function
(xTree2, # has cpName(choppedDat), varGenes
 washedDat, # has expDat
 stDat,
 propTrain=0.25,
 nTrees=200,
 numGenes = 500) #number of genes for building the classifiers
{
  #can't get out the nested loop structure 
  
  #check to see how many levels there are
  tmp <- as.list(xTree2) #convert xTree2 into tmp
  classifiers_list <- list()
  geneLists <- list()
  
  #top level
  topL <- xTree2$grp_list[[2]]
  stDat_topL <- cbind(stDat, topL)
  diffExp <- gnrAll(washedDat[['expDat']], stDat_topL$topL)
  cts<-as.vector(unique(stDat_topL[,"topL"]))
  
  for(i in 1:length(cts)){
    predictors <- getClassGenes(diffExp[[i]],numGenes)
    geneLists[[cts[i]]] <- predictors
  }
  # Make classifiers
  cat("Making classifiers (this can take awhile) ...\n")
  
  myRFs<-makeRFs(washedDat[['expDat']], stDat_topL, geneLists, dLevel="topL", nTrees=nTrees)

  classifiers_list[["topLevel"]] <- list(classResVal=list(classRes=sc_classify(myRFs, washedDat[['expDat']], geneLists), stTrain=stDat_topL),
                                         classifiers=myRFs,
                                         predictors=geneLists)
  
  #bottome lvel
  if(length(xTree2$grp_list) > 2){
    bottomL <- xTree2$grp_list[[length(xTree2$grp_list)]]
    stDat_bottomL <- cbind(stDat, bottomL)
    diffExp <- gnrAll(washedDat[['expDat']], stDat_bottomL$bottomL)
    cts<-as.vector(unique(stDat_bottomL[,"bottomL"]))
    
    for(i in 1:length(cts)){
      predictors <- getTopGenes(diffExp[[i]],numGenes)
      geneLists[[cts[i]]] <- predictors
    }
    
    # Make classifiers
    cat("Making classifiers (this can take awhile) ...\n")
    
    myRFs<-makeRFs(washedDat[['expDat']], stDat_bottomL, geneLists, dLevel="bottomL",
                   nTrees=nTrees)

    classifiers_list[["bottomLevel"]] <- list(classResVal=list(classRes=sc_classify(myRFs, washedDat[['expDat']], geneLists),
                                                               stTrain=stDat_bottomL),
                                              classifiers=myRFs,
                                              predictors=geneLists)
    
  }     
  
  return(classifiers_list)  
}

#' getClassGenes
#'
#' extract genes for training classifier
#' @param diffRes a df with cval, holm, rownames=genes
#' @param topX number of genes to select
#' @param bottom boolean if ture use the top x genes with - cvals
#'
#' @return vector of genes
#'
#' @export
getClassGenes<-function(
 diffRes,
 topX=25,
 bottom=TRUE)
 {
   #exclude NAs
   xi<-which(!is.na(diffRes$cval))
   diffRes<-diffRes[xi,]   
   diffRes<-diffRes[order(diffRes$cval, decreasing=TRUE),]
   ans<-rownames(diffRes[1:topX,])
   if(bottom){
     ans<-append(ans, rownames( diffRes[nrow(diffRes) - ((topX-1):0),]))
   }
   ans
 }
