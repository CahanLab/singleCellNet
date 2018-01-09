# singleCellNet
# (C) Patrick Cahan 2012-2017

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




