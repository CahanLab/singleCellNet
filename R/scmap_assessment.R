#build assessment pipeline for scmap

assessmentReport_scmap <- function(siml, #similarity from scmap
                                  true_label, #sample table where cells in query are in the training with name 
                                  pred_label,
                                  resolution = 0.005,# increment at which to evalutate classification
                                  classLevels = "description2",
                                  dLevelSID = "sample_name"){
  
  report <- list()
  
  #cohen's kappa, accuracy
  
  cm = as.matrix(table(Actual = true_label, Predicted = pred_label))
  
  #in case of misclassfication where there are classifiers that are not used
  if(length(setdiff(unique(true_label), unique(pred_label))) != 0){
    misCol <- setdiff(unique(true_label), unique(pred_label))
    for(i in 1:length(misCol)){
      cm <- cbind(cm, rep(0, nrow(cm)))
    }
    colnames(cm)[(ncol(cm) - length(misCol) +1) : ncol(cm)] <- misCol
  }
  
  if(length(setdiff(unique(pred_label), unique(true_label))) != 0){
    misRow <- setdiff(unique(pred_label), unique(true_label))
    for(i in 1:length(misRow)){
      cm <- rbind(cm, rep(0, ncol(cm)))
    }
    rownames(cm)[(nrow(cm) - length(misRow) +1) : nrow(cm)] <- misRow
  }
  
  cm <- cm[,colnames(cm)[match(rownames(cm),colnames(cm))]]  
  report[['cm']] <- cm
  
  #sort table names accordigly
  
  n = sum(cm) # number of instances
  nc = nrow(cm) # number of classes
  diag = diag(cm) # number of correctly classified instances per class 
  rowsums = apply(cm, 1, sum) # number of instances per class
  colsums = apply(cm, 2, sum) # number of predictions per class 
  p = rowsums / n # distribution of instances over the actual classes
  q = colsums / n # distribution of instances over the predicted classes
  expAccuracy = sum(p*q)
  accuracy = sum(diag) / n 
  report[['kappa']] <- (accuracy - expAccuracy) / (1 - expAccuracy)
  report[['accuracy']] <- accuracy
  
  return(report)
  
}

plot_scmapAssess <- function(assessed, method = "scmap"){
 metric <- matrix(0, ncol = 2, nrow = 1)
  colnames(metric) <- c("cohen's kappa", "accuracy")
  rownames(metric) <- "value"
  metric[,1:2] <- c(assessed$kappa, assessed$accuracy) 
  metric <- as.data.frame(metric)
  
  p1<-ggplot(metric, aes(x="cohen's kappa", y = metric[1,1])) + geom_bar(stat="identity") +xlab("") + ylab("") + theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +  ylim(0,1) + theme(legend.position="none")
  
  p2<-ggplot(metric, aes(x="accuracy", y = metric[1,2])) + geom_bar(stat="identity") +xlab("") + ylab("") + theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) + ylim(0,1) + theme(legend.position="none")
  
  p1 + p2
  
}



prep_SimilarityMatrix <- function(projection = sce_test,
                                  index_list =  metadata(sce)$scmap_cluster_index, #please compare one index_list at a time
                                  expTest = expTest #only used for the sample names
                                  ){
  
  answer <- list()
  labels <- list()
  simls <- list()
  
  index <- index_list
  # find and select only common features, then subset both datasets
  tmp <- setFeatures(projection, rownames(index))
  index <- index[rownames(index) %in% rowData(tmp)$feature_symbol[rowData(tmp)$scmap_features], , drop = FALSE]
  tmp <- tmp[rowData(tmp)$scmap_features, ]
  
  if (nrow(index) < 10) {
    warning("There are less than ten features in common between the `reference` and `projection` datasets. Most probably they come from different organisms! Please redefine your query!")
    return(projection)
  }  
  
  # get expression values of the projection dataset
  proj_exprs <- as.matrix(logcounts(tmp))
  rownames(proj_exprs) <- rowData(tmp)$feature_symbol
  
  # prepare projection dataset
  proj_exprs <- proj_exprs[order(rownames(proj_exprs)), ]
  
  # calculate similarities and correlations
  tmp <- t(index)
  res_cosine <- proxy::simil(tmp, t(proj_exprs), method = "cosine")
  res_cosine <- matrix(res_cosine, ncol = nrow(tmp), byrow = TRUE) #cols are training cell types; rows are cells
  colnames(res_cosine) <- colnames(index)
  rownames(res_cosine) <- colnames(expTest)
  histogram(rowSums(res_cosine)) #show distribution of rowSum for res_cosine
  res_cosine <- res_cosine/rowSums(res_cosine) 
  
  res_pearson <- cor(index, proj_exprs, method = "pearson")
  res_pearson <- matrix(res_pearson, ncol = nrow(tmp), byrow = TRUE)  #cols are training cell types; rows are cells
  colnames(res_pearson) <- colnames(index)
  rownames(res_pearson) <- colnames(expTest)
  histogram(rowSums(res_pearson)) #show distribution of rowSum for res_cosine
  res_pearson <- res_pearson/rowSums(res_pearson)
  
  res_spearman <- cor(index, proj_exprs, method = "spearman")
  res_spearman <- matrix(res_spearman, ncol = nrow(tmp), byrow = TRUE)  #cols are training cell types; rows are cells
  colnames(res_spearman) <- colnames(index)
  rownames(res_spearman) <- colnames(expTest)
  histogram(rowSums(res_spearman))
  res_spearman <- res_spearman/rowSums(res_spearman)
  
  return(answer= list(res_cosine = res_cosine, res_pearson = res_pearson, res_spearman = res_spearman))
}
