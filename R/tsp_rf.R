# for TSP-RF
# (C) 2018 Patrick Cahan

#' find best pairs
#'
#' find best pairs
#'
#' @param expDat expDat
#' @param cellLabels named vector of cell groups
#'
#' @return vector of pairs
#' 
#' @export
gnrBP<-function(
  expDat,
  cellLabels,
  topX=50){

	ans<-vector()
    myPatternG<-sc_sampR_to_pattern(as.character(cellLabels))
    for(i in seq(length(myPatternG))){
    	cat(i,"\n")
    	xres<-sc_testPattern(myPatternG[[i]], expDat=expDat)
    	tmpAns<-findBestPairs(xres, topX)
    	ans<-append(ans, tmpAns)
    }
    unique(ans)
}


#' find candidate classifier-worthy genes
#'
#' find candidate classifier-worthy genes
#'
#' @param expDat expDat
#' @param sampTab sampTab
#' @param dLevel dLevel
#' @param topX topX
#' @param dThresh dThresh
#' @param alpha1 alpha1
#' @param alpha2 alpha2
#' @param mu mu
#'
#' @return list of cgenes and grps
#' 
#' @export
findClassyGenes<-function
(expDat,
 sampTab,
 dLevel,
 topX=25,
 dThresh=0,
 alpha1=0.05,
 alpha2=.001,
 mu=2)
{
	gsTrain<-sc_statTab(expDat, dThresh=dThresh)
	ggenes<-sc_filterGenes(gsTrain, alpha1=alpha1, alpha2=alpha2, mu=mu)
	grps<-as.vector(sampTab[,dLevel])
	names(grps)<-rownames(sampTab)
	xdiff<-gnrAll(expDat[ggenes,], grps)
	cgenes<-lapply(xdiff, getClassGenes, topX=topX)
	cgenes<-unique(unlist(cgenes))
	list(cgenes=cgenes, grps=grps)
}


#hm_gpa_sel(expTrain, cgenes, grps, maxPerGrp=25, toScale=T, cRow=F, cCol=F,font=4)


#' makes complete gene-to-gene comparison
#'
#' @param expDat expDat
#'
#' @return list of list(genes=data frame, expDat=binary matrix)
#'
#' @export
pair_transform<-function # convert to a vector of length = length(vect)^2 - 1 /2
(expDat){
	ngenes<-nrow(expDat)
	genes<-rownames(expDat)
	ans<-matrix(0, nrow=ngenes*(ngenes-1)/2, ncol=ncol(expDat))
	pair_index<-1
	genes1<-vector()
	genes2<-vector()
	for(i in 1:ngenes){
		for(j in 1:ngenes){
			if(j>i){
				genes1<-append(genes1, genes[i])
				genes2<-append(genes2, genes[j])
				ans[pair_index,]<-as.numeric(expDat[i,]>expDat[j,]) 
				pair_index<-pair_index +1
			}
		}
	}
	colnames(ans)<-colnames(expDat)
	tList2 <- list(genes=data.frame(g1=genes1, g2=genes2), tDat=ans)

	pairNames<-paste(tList2[[1]][,1], "_",tList2[[1]][,2], sep='')

	pairDat<-tList2[[2]]
	rownames(pairDat)<-pairNames
	pairDat
}


#' makes complete gene-to-gene comparison
#'
#' @param expDat expDat
#' @param genePairs genePairs
#'
#' @return matrix indicating which gene of a pair is greater
#'
#' @export
query_transform<-function # convert to a vector of length = length(vect)^2 - 1 /2
(expDat,
 genePairs #vector of strings indicating pairs to compare
 ){
	genes<-strsplit(genePairs, "_")
    ans<-matrix(0, nrow=length(genes), ncol=ncol(expDat))
	pair_index<-1
	genes1<-vector()
	genes2<-vector()
	for(i in seq(length(genes))){
		ans[i,]<-as.numeric(expDat[genes[[i]][1],]>expDat[genes[[i]][2],])
	}
	colnames(ans)<-colnames(expDat)
	rownames(ans)<-genePairs
	ans
}

#' finds the best pairs to use
#'
#' @param xdiff xdiff
#' @param n number of pairs
#' @param maxPer indicates the number of pairs that a gene is allowed to be in
#'
#' @return vector of good pairs
#'
#' @export
findBestPairs<-function # find best and diverse set of pairs
(xdiff,
 n=50,
 maxPer=3){
 
 	
 	xdiff<-xdiff[order(xdiff$cval, decreasing=TRUE),]
 	genes<-unique(unlist(strsplit(rownames(xdiff), "_")))
 	countList <- rep(0, length(genes))
 	names(countList) <- genes	

 	i<-0
 	ans<-vector()
 	xdiff_index<-1
 	pair_names<-rownames(xdiff)
 	while(i < n ){
 		tmpAns<-pair_names[xdiff_index]
 		tgp <- unlist(strsplit(tmpAns, "_"))
 		if( (countList[ tgp[1] ] < maxPer) & (countList[ tgp[2] ] < maxPer )){
 			ans<-append(ans, tmpAns)
 			countList[ tgp[1] ] <- countList[ tgp[1] ]+ 1
 			countList[ tgp[2] ] <- countList[ tgp[2] ]+ 1
 			i<-i+1
 		}
 		xdiff_index <- xdiff_index + 1
 	}
 	ans
 }


#' @export
addRandToSampTab<-function(classRes, sampTab, desc, id="cell_name"){
	cNames<-colnames(classRes)
	snames<-rownames(sampTab)

	rnames<-setdiff(cNames, snames)
	cat("number of random samples: ",length(rnames), "\n")

	stNew<-data.frame(rid=rnames, rdesc=rep("rand", length(rnames)))
	stTop<-sampTab[,c(id, desc)]
	colnames(stNew)<-c(id, desc)
	ans<-rbind(stTop, stNew)
	rownames(ans)<-colnames(classRes)
	ans
}














