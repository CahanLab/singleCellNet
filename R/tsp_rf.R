# for TSP-RF
# (C) 2018 Patrick Cahan



#' enable cross-species comparison
#'
#' rename and subset query expDat to orthologs in training data
#'
#' @param expQuery expQuery
#' @param expTrain expTrain
#' @param orthTable orthTable
#' @param speciesQuery speciesQuery default human
#' @param speciesTrain speciesTrain default mouse
#' 
#' @return list (expQuery, expTrain)
#' 
#' @export
csRenameOrth<-function(
	expQuery,
	expTrain,
	orthTable,
	speciesQuery='human',
	speciesTrain='mouse'){

 	# what genes in the orth table are in the query table?
	rownames(orthTable)<-as.vector(orthTable[,speciesQuery])
	cgenes<-intersect(rownames(expQuery), rownames(orthTable))
	cat("query genes in ortholog table = ",length(cgenes),"\n")

	# of these, which are in training data?
	oTab<- orthTable[cgenes,]
	rownames(oTab) <- as.vector(oTab[,speciesTrain])
	ccGenes<- intersect(rownames(expTrain), rownames(oTab))
	cat("training genes in ortholog table and query data = ",length(ccGenes),"\n")
	# should put a check here for sufficient number of genes 
	oTab<-oTab[ccGenes,]
	qGenes<-as.vector(oTab[,speciesQuery])
	expQuery <- expQuery[qGenes,]
	rownames(expQuery) <- rownames(oTab) 
	expTrain <- expTrain[rownames(oTab),]
	list(expQuery=expQuery, expTrain=expTrain)
}


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
	cgenes2<-unique(unlist(cgenes))
	list(cgenes=cgenes2, grps=grps, by_ct = cgenes)
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

#' find genes higher in a cluster compared to all other cells
#'
#' ind genes higher in a cluster compared to all other cells
#'
#' @param expDat expDat
#' @param cellLabels named vector of cell groups
#'
#' @return list of diffExp data framnes
#' 
#' @export
gnrAll<-function(
  expDat,
  cellLabels){

  myPatternG<-sc_sampR_to_pattern(as.character(cellLabels))
  # sparse matrix?
  if(class(expDat)[1]!='matrix'){
    expTrans = Matrix::t(expDat)
  }
  else{
    expTrans = t(expDat)
  }
  specificSets<-lapply(myPatternG, sc_testPatternTrans, expDat=expTrans)
  cat("Done testing\n")

#  grpOrder<-myGrpSort(cellLabels)

#  specificSets[grpOrder]

  specificSets
}


makePairTab<-function(genes){
	pTab<-t(combn(genes, 2))
	colnames(pTab)<-c("genes1", "genes2")
	pTab<-cbind(pTab, pairName=paste(pTab[,1], "_",pTab[,2], sep=''))
	pTab
}



#' makes vector of gene pairs, iterates over this and computes pairDat, sc_testPattern, then, at the end, findBestPairs 
#'
#' @param expDat expDat
#' @param cell_labels named vector, value is grp, name is cell name
#' @param topX 50
#'
#' @return vector of gene-pair names
#'
#' @export
ptGetTop<-function
(expDat,
 cell_labels,
 topX=50,
 sliceSize = 5e3,
 ncores = detectCores()){

	ans<-vector()
	genes<-rownames(expDat)

	mcCores <- 1
	if(ncores>1){
		mcCores <- ncores - 1
	}
	cat(ncores, "  --> ", mcCores,"\n")

	# make a data frame of pairs of genes that will be sliced later
	if(FALSE){
		cat("Making pairTable\n")
		ngenes<-nrow(expDat)
		genes<-rownames(expDat)
		genes1<-vector()
		genes2<-vector()
		for(i in 1:ngenes){ # replace with combn?
			for(j in 1:ngenes){
				if(j>i){
					genes1<-append(genes1, genes[i])
					genes2<-append(genes2, genes[j])				
				}
			}
		}

		pairTab = data.frame(genes1=genes1, genes2=genes2)
		pairNames<-paste(pairTab[,1], "_",pairTab[,2], sep='')
		pairTab <- cbind(pairTab, pairName=pairNames)
	}

	pairTab<-makePairTab(genes)

	###
	# setup tmp ans list of sc_testPattern
	cat("setup ans and make pattern\n")
	grps<-unique(cell_labels)
	myPatternG<-sc_sampR_to_pattern(as.character(cell_labels))
	statList<-list()
	for(grp in grps){
		statList[[grp]]<-data.frame()
	}

	# make the pairedDat, and run sc_testPattern
	cat("make pairDat on slice and test\n")
	nPairs = nrow(pairTab)
	cat("nPairs = ",nPairs,"\n")
	str = 1
	stp = min(c(sliceSize, nPairs))
	while(str <= nPairs){
		if(stp>nPairs){
			stp <- nPairs
		}
		cat(str,"-", stp,"\n")
		tmpTab<-pairTab[str:stp,]
		tmpPdat<-ptSmall(expDat, tmpTab)

	### new
	
	if (Sys.info()[['sysname']] == "Windows") {
   		 tmpAns<-lapply(myPatternG, sc_testPattern, expDat=tmpPdat)
  	}
  	else {
    		tmpAns<-parallel::mclapply(myPatternG, sc_testPattern, expDat=tmpPdat, mc.cores=mcCores) # this code cannot run on windows
  	}
	

	for(gi in seq(length(myPatternG))){
	    	grp<-grps[[gi]]
	    	#cat(i, " grp: ",grp,"\n")
    		statList[[grp]]<-rbind( statList[[grp]],  tmpAns[[grp]])
    	}


    	str = stp+1
   	 	stp = str + sliceSize - 1		
	}

	cat("compile results\n")
	for(grp in grps){
    	tmpAns<-findBestPairs(statList[[grp]], topX)
    	ans<-append(ans, tmpAns)
    }
    unique(ans)
}

ptSmall<-function
(expDat,
 pTab){
	npairs = nrow(pTab)
	ans<-matrix(0, nrow=npairs, ncol=ncol(expDat))
	genes1<-as.vector(pTab[,"genes1"])
	genes2<-as.vector(pTab[,"genes2"])

    for(i in seq(nrow(pTab))){
    	#cat(genes1[i], ": ", genes2[i],"\n")
    	ans[i,]<-as.numeric(expDat[genes1[i],]>expDat[genes2[i],]) 
    }
	colnames(ans)<-colnames(expDat)
	rownames(ans)<-as.vector(pTab[,"pairName"])
	ans
}


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














