
library(CellNet)
library(singleCellNet)
library(SC3)
library(scater)
library(Rtsne)
library(pheatmap)
library(plyr)
#this version sepeartes group assignment in a sepearte name vector and replace that in sampTab before it passes in the pipe_butter function

#steamed
steam_sc3 <- function(
  dataMat,
  sampTab,
  k_estimator = FALSE, #optimal k
  ks = NULL,  #k mean range, it will be messed up if k starts with 1
  silhouette_threshold = 0.8,
  pct_dropout_min = 10, 
  pct_dropout_max = 90,
  gene_filter = FALSE,
  d_region_min = 0.04, #range of the epigenevector 
  d_region_max = 0.07, 
  svm_num_cells = NULL, 
  svm_train_inds = NULL, 
  svm_max = 5000, 
  n_cores = NULL, 
  kmeans_nstart = NULL, 
  kmeans_iter_max = 1e+09,  
  biology = FALSE, 
  rand_seed = 1
){
  
  #convert expression matrix into in sceset object
  object <- newSCESet(exprsData = dataMat)
  ans<-NULL
  if(k_estimator){
    cat("calculating the optimal k value\n")
    object <- sc3_estimate_k(object)
    optimal_k <- object@sc3$k_estimation
    
    cat("clustering\n")
    object <- sc3(object, ks = optimal_k, gene_filter = FALSE)
    
    #return the clusting assignment, the column of sc3_'k'_cluster contains the annotation information
    group_list <- object@phenoData@data
    group_list <- group_list[match(rownames(sampTab), rownames(group_list)),]
    rownames(group_list) <- rownmaes(sampTab)
    colnames(group_list) <- as.character(optimal_k)
    
    #convert cluster assignment to characters
    group_list <- data.frame(lapply(group_list, as.character), stringsAsFactors=FALSE)
    
    opt_params <- optimal_k
    args <- as.list(match.call())
    ans<-list(sampTab = sampTab, args = args, opt_params = opt_params, group_list = group_list)
    
  } else {
    cat("sc3 clustering\n")
    opt_params <- rep(NA, 20) 
    group_list <- data.frame(matrix(nrow = nrow(sampTab), ncol =0))
    
    #calculate average silhouette index, set the threshold to be 0.7 
    #collect the ks that passes that threshold 
    for (k in ks) {
      object <- sc3(object, ks = k, gene_filter = FALSE)
      index <- object@sc3$consensus[[1]]
      index_sum <- summary(index$silhouette, FUN = mean)
      index_avg <- sum(index_sum[[2]])/k
      if (index_avg >= silhouette_threshold){
        #get the group assignment
        tmp_ans <- object@phenoData@data
        colnames(tmp_ans) <- as.character(k)
        group_list <- cbind(group_list, tmp_ans)
        opt_params[k] <- k
      } 
    }
    
    #return the clusting assignment, the column of sc3_'k'_cluster contains the annotation information
    opt_params <- opt_params[!is.na(opt_params)]
    
    if(!empty(group_list)){ #check to see with the situation where group_list is empty
      #need to order the grouplist according to the sampTab here
      group_list$names<-rownames(group_list)
      group_list <- group_list[match(rownames(sampTab), rownames(group_list)),]
      
      #convert cluster assignment to characters
      group_list <- data.frame(lapply(group_list, as.character), stringsAsFactors=FALSE)
      group_list<-as.data.frame(group_list[,-ncol(group_list)])
      
      rownames(group_list) <- rownames(sampTab)
      colnames(group_list) <- opt_params
    }
    
    if(empty(group_list)){
      cat("group_list is empty, you may want to re-adjust silhouette_threshold\n")
    }
    
    args <- as.list(match.call())
    ans<-list(sampTab = sampTab, args = args, opt_params = opt_params, group_list = group_list)
  }
  
  return(ans)
}

#to make pipeSteamed object
pipe_sc3 <- function
(washedDat,
 sampTab,
 silhouette_threshold = 0.8,
 ks=NULL,
 topPC=20,
 zThresh = 2)
{
  cp_pca<-chop_pca(washedDat[['expDat']], washedDat[['geneStats']], 
                   zThresh=zThresh, meanType="overall_mean")
  
  # tsne
  cat("tsne-ing\n")
  cp_tsne<-chop_tsne(cp_pca[['choppedDat']][,1:topPC], perplexity=30, theta=.3)
  
  #sc3
  cat("SC3 clustering\n")
  stm_sc3 <- steam_sc3(washedDat[['expDat']], sampTab, k_estimator = FALSE, ks = ks, gene_filter = FALSE,silhouette_threshold = silhouette_threshold, biology = FALSE)
  
  if(empty(group_list)) {
    cat("group_list is empty\n")
    cat("please adjust the silhouette threshold and remake your pipeSteamed object\n")
    
  }
  
  #there is group_list in steam_sc3
  list(cp_pca = cp_pca, cp_tsne = cp_tsne, steamed = stm_sc3)
}


#access the classifer performance
pipe_cAss_sc3 <- function
(washedDat,
 sampTab,
 ks = NULL,
 silhouette_threshold = 0.8,
 topPC = 20,
 zThresh = 2)
{
  # all methods need PCA dimension reduction, this will gives us varGenes
  cat("reducing dimensionality\n")
  cp_pca<-chop_pca(washedDat[['expDat']], washedDat[['geneStats']], 
                   zThresh=zThresh, meanType="overall_mean")
  
  # tsne
  cat("tsne-ing\n")
  cp_tsne<-chop_tsne(cp_pca[['choppedDat']][,1:topPC], perplexity=30, theta=.3)
  
  #sc3
  cat("SC3 clustering\n")
  stm_sc3 <- steam_sc3(washedDat[['expDat']], sampTab, k_estimator = FALSE, ks = ks, silhouette_threshold = silhouette_threshold, gene_filter = FALSE, biology = FALSE)
  
  #pre-butter each, to split and divide up the training data
  cat("pre-butter sc3\n")
  #this is where we start to split the group list and adding it to the sampTab
  stDat <- stm_sc3[['sampTab']]
  group_list <- stm_sc3[['group_list']]
  
  if(empty(group_list)) {
    cat("group_list is empty\n")
    stop("please adjust the silhouette threshold\n")   
  }
  
  opt_params <- stm_sc3[['opt_params']]
  divide_sampTab_list <- list()
  
  #instead of straight appending, just write a for loop to add the assignment
  for (i in 1:ncol(group_list)) {
    k = colnames(group_list)[i]
    divide_sampTab_name <- paste0("classRes_sc3_", k, sep = "")
    stDat_tmp <- cbind(sampTab, group_list[,i])
    colnames(stDat_tmp)[ncol(stDat_tmp)] <- "group"
    divide_sampTab_list[[divide_sampTab_name]]<-prebutter(cp_pca$varGenes, washedDat[['expDat']], stDat_tmp,
                                                          propTrain=.25, partOn="group", nTrees=200)
    
  }
  
  #assess and iterate through the classRes_sc3 list, name them accordingly
  #for loop to iterate through the classifier_list
  classRes_sc3_list <- list()
  
  for (i in 1:length(opt_params)) {
    k = opt_params[i]
    classRes_name <- paste0("cAss_sc3_", k, sep = "")
    
    #iterate through the divide_sampTab_list to identify the classRes_sc3 per k
    divide_sampTab_name <- paste0("classRes_sc3_", k, sep = "")
    classRes_sc3<- divide_sampTab_list[[divide_sampTab_name]]
    
    #classifying/access and put into a list
    classRes_sc3_list[[classRes_name]] <- easy_assess(classRes_sc3[['classRes']], classRes_sc3[['stVal']])
    
  }
  
  return(classRes_sc3_list)
}


#make classifiers
pipe_butter_sc3<-function
(pipeSteamed, # has cpName(choppedDat), varGenes
 washedDat, # has expDat
 propTrain=0.25,
 partOn="group",
 nTrees=200)
{
  geneLists<-list()
  #this is where i need to loop throught the divide sampTab loop
  stDat<-pipeSteamed[['steamed']][['sampTab']]
  group_list <- pipeSteamed[['steamed']][['group_list']]
  
  if(empty(group_list)) {
    cat("group_list is empty\n")
    stop("please adjust the silhouette threshold and remake your pipeSteamed object\n")   
  }
  
  opt_params <- pipeSteamed[['steamed']][['opt_params']]
  predictors<-pipeSteamed[['cp_pca']][['varGenes']]
  classResVal_list <- list()
  
  for( i in 1: length(opt_params)){
    stDat_tmp <- cbind(stDat, group_list[,i])
    colnames(stDat_tmp)[ncol(stDat_tmp)] <- "group"
    k = opt_params[i]
    
    cts<-as.vector(unique(stDat_tmp[,partOn]))
    for(ct in cts){
      geneLists[[ct]]<-predictors
    }
    
    # split into training and test data
    ttList<-divide_sampTab(stDat_tmp, propTrain, partOn)
    stTrain<-ttList[['stTrain']]
    
    # make RFs
    myRFs<-makeRFs(washedDat[['expDat']][predictors,rownames(stTrain)], stTrain, geneLists, dLevel=partOn,
                   nTrees=nTrees)
    
    # classify held out data
    stVal<-ttList[['stVal']]
    
    classResVal_name <- paste0("classResVal_", k, sep = "")
    
    classResVal_list[[classResVal_name]] <- list(classResVal=list(classRes=sc_classify(myRFs, washedDat[['expDat']][predictors,rownames(stVal)], geneLists),
                                                                  stVal=stVal),
                                                 classifiers=myRFs,
                                                 predictors=geneLists)
  }
  
  return(classResVal_list)
}



