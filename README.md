# singleCellNet

### Introduction
See [CellNet](https://github.com/pcahan1/CellNet) for an introduction to CellNet, how to use it on bulk RNA-Seq, and how to analyze single cell RNA-Seq (scRNA-Seq) data with classifiers trained on bulk RNA-Seq. Here, we illustrate

- how to build and assess single cell classifiers

- how to use these classifiers to quantify 'cell identity' from query scRNA-Seq data

- how to cluster scRNA-Seq data using our 'cluster by competition' method

### DATA

In this example, we use a subset of the Tabula Muris data to train singleCellNet. To learn more about the Tabula Muris project, see the [manuscript])(https://www.biorxiv.org/content/early/2018/03/29/237446). As query data, we use scRNA-Seq of kidney cells as reported in [Park et al 2018](https://www.ncbi.nlm.nih.gov/pubmed/29622724). You can download this data here:

| APPLICATION | METADATA | EXPRESSION |
|-------------|----------|------------|
| Query       | [metadata](https://s3.amazonaws.com/cnobjects/singleCellNet/examples/sampTab_Park_MouseKidney_062118.rda) | [expression data](https://s3.amazonaws.com/cnobjects/singleCellNet/examples/GSE107585_Mouse_kidney_single_cell_datamatrix.txt.gz) |
| Training    | [metadata](https://s3.amazonaws.com/cnobjects/singleCellNet/examples/sampTab_TM_053018.rda) | [expression data](https://s3.amazonaws.com/cnobjects/singleCellNet/examples/expTM_Raw_053018.rda) |

N.B. The query expression data needs to be decompressed before loading it into R. 


#### Setup
```R
library(fgsea)
library(devtools)
install_github("pcahan1/singleCellNet", ref="tsp_rf_pc", auth="your_token")
library(singleCellNet)

library(RColorBrewer)
library(pheatmap)
library(randomForest)
library(viridis)
library(ggplot2)

mydate<-utils_myDate()
```

#### Fetch the data (optional if you have alread done this)
```R
download.file("https://s3.amazonaws.com/cnobjects/singleCellNet/examples/sampTab_Park_MouseKidney_062118.rda", "sampTab_Park_MouseKidney_062118.rda")

download.file("https://s3.amazonaws.com/cnobjects/singleCellNet/examples/expDat_Park_MouseKidney_062218.rda", "expDat_Park_MouseKidney_062218.rda")

download.file("https://s3.amazonaws.com/cnobjects/singleCellNet/examples/expTM_Raw_053018.rda", "expTM_Raw_053018.rda")

download.file("https://s3.amazonaws.com/cnobjects/singleCellNet/examples/sampTab_TM_053018.rda", "sampTab_TM_053018.rda")
```

#### Load query data
```R
stPark<-utils_loadObject("sampTab_Park_MouseKidney_062118.rda")
expPark<-utils_loadObject("expDat_Park_MouseKidney_062218.rda")
dim(expPark)
[1] 16272 43745

genesPark<-rownames(expPark)
```

#### Load the training data
```R
expTMraw<-utils_loadObject("expTM_Raw_053018.rda")
dim(expTMraw)
[1] 23433 24936

stTM<-utils_loadObject("sampTab_TM_053018.rda")
dim(stTM)
[1] 24936    17

stTM<-droplevels(stTM)
```

#### Find genes in common to the data sets
```R
commonGenes<-intersect(rownames(expTMraw), genesPark)
length(commonGenes)
[1] 13831
```


#### Normalize the training data
```R
expTMnorm<-trans_prop(weighted_down(expTMraw[commonGenes,], 1.5e3), 1e4)
```

#### Find the best set of classifier genes
```R
stList<-splitCommon(stTM, ncells=100, dLevel="newAnn")
stTrain<-stList[[1]]
expTrain<-expTMnorm[,rownames(stTrain)]

system.time(cgenes2<-findClassyGenes(expTrain, stTrain, "newAnn", topX=10))
   user  system elapsed 
 92.374  19.027 111.651 

cgenesA<-cgenes2[['cgenes']]
grps<-cgenes2[['grps']]
length(cgenesA)
[1] 481

# heatmap these genes
hm_gpa_sel(expTrain, cgenesA, grps, maxPerGrp=5, toScale=T, cRow=F, cCol=F,font=4)
```
<img src="md_img/hm_tabulaMuris.png">

#### TSP transform the training data
```R
system.time(pairDat<-pair_transform(expTrain[cgenesA,]))
   user  system elapsed 
 83.668  30.796 114.485
 ```

####  Find the best pairs
 ```R
 system.time(xpairs<-gnrBP(pairDat, grps)) # will take ~ 35 minutes
    user   system  elapsed 
1662.301  413.081 2124.899 

length(xpairs)
[1] 1562
```

#### Train the classifier
```R
system.time(rf_tspAll<-sc_makeClassifier(pairDat[xpairs,], genes=xpairs, groups=grps, nRand=100, ntrees=1000))    
   user  system elapsed 
700.098   1.090 701.430
```


#### Apply to held out data -- this is the place to add the multi-class assessment
```R
stTest<-stList[[2]]

system.time(expQtransAll<-query_transform(expTMraw[cgenesA,rownames(stTest)], xpairs))
 
    user  system elapsed 
  4.221   2.751  11.369 


nrand<-100
system.time(classRes_val_all<-rf_classPredict(rf_tspAll, expQtransAll, numRand=nrand))
  user  system elapsed 
 37.136   1.522  38.691 

sla<-as.vector(stTest$newAnn)
names(sla)<-rownames(stTest)
slaRand<-rep("rand", nrand)
names(slaRand)<-paste("rand_", 1:nrand, sep='')

sla<-append(sla, slaRand)

# heatmap classification result
sc_hmClass(classRes_val_all, sla, max=300, isBig=TRUE)
```
<img src="md_img/hmClass_validation.png">

#### Apply to Park et al query data
```R
system.time(kidTransAll<-query_transform(expPark[cgenesA,], xpairs))
   user  system elapsed 
  8.594   3.314  15.704 
  
nqRand<-100

system.time(crParkall<-rf_classPredict(rf_tspAll, kidTransAll, numRand=nqRand))
  user  system elapsed 
 78.520   3.873  82.49

sgrp<-as.vector(stPark$description1)
names(sgrp)<-rownames(stPark)
grpRand<-rep("rand", nqRand)
names(grpRand)<-paste("rand_", 1:nqRand, sep='')
sgrp<-append(sgrp, grpRand)

# heatmap classification result
sc_hmClass(crParkall, sgrp, max=5000, isBig=TRUE, cCol=F, font=8)
```
<img src="md_img/hmClass_Park.png">


### Skyline plot of classification results
```R
stKid2<-addRandToSampTab(crParkall, stPark, "description1", "sample_name")
skylineClass(crParkall, "T cell", stKid2, "description1",.25, "sample_name")
```
<img src="md_img/skyline_Tcell_Park.png">



#### Determine cell groups yourself using Cluster by Competition (CBC)
```R
library(cluster)
library(pcaMethods)
library(rpca)
library(data.tree)
library(dbscan)
library(Rtsne)

gstats<-sc_statTab(expTrain, dThresh=1)
ggenes<-sc_filterGenes(gstats, alpha1=0.025, alpha2=.001, mu=2)
length(ggenes)
[1] 7777

system.time(xTree_auto<-gpa_recurse(expTrain[ggenes,], zThresh=2, maxLevel=3, nPCs=2, SilhDrop=0.5, methods=c("cutree","sNN_clust","kmeans"), dThresh=1, pcaMethod="prcomp",k=2:8, minClusterSize=20, silMin=FALSE, pcAuto=TRUE))

   user  system elapsed 
 77.129   4.513  81.654 

system.time(ts3<-pca_to_tsne(expTrain[ggenes,], xTree_auto, perplexity=30, theta=0.25, weighted=FALSE))
 user  system elapsed 
 44.945   0.410  45.370


stT1<-stTrain
stT1<-cbind(stT1, cluster=xTree_auto$groups)
sampTab<-cbind(sampTab, cluster2=xTree_auto$grp_list[[2]])
sampTab<-cbind(sampTab, cluster3=xTree_auto$grp_list[[3]])

plot_tsne(stT1, ts3, cname="cluster")
```
<img src="md_img/tsne_TM_cluster.png">

#### Black background; good for presentations
```R
plot_tsne(stT1, ts3, cname="cluster", themeWhite=FALSE)
```
<img src="md_img/tsne_TM_cluster_black.png">


Color by category provided by Tabula Muris
```R
plot_tsne(stT1, ts3, cname="newAnn")
```
<img src="md_img/tsne_TM_label.png">

 
#### Differential expression
```R
system.time(xdiff<-gnrAll(expTrain[ggenes,], xTree_auto$groups))
   user  system elapsed 
 71.072  16.591  87.684 

x1<-lapply( xdiff, getTopGenes, 5)
newOrder<-reorderCellsByGrp( xTree_auto$groups , names(x1))

hm_gpa_sel(expTrain, c(unique(unlist(x1))), newOrder, maxPerGrp=20, toScale=T, cRow=F, cCol=F,font=5)
```
<img src="md_img/hm_diffExp_TM.png">

#### tsne plot specific genes 
```R 
tsneMultsimp(ts3, expTrain, c("Ear2","Tnnt2","Cd3e","Apoc3"))
```
<img src="md_img/tsne_4_genes.png">

#### tsne plot specific genes; black background
```R 
tsneMultsimp(ts3, expTrain, c("Ear2","Tnnt2","Cd3e","Apoc3"), revCol=FALSE,themeWhite=FALSE)
```
<img src="md_img/tsne_4_genes_black.png">

#### Gene set enrichment analysis

Download and install fgsea (https://github.com/ctlab/fgsea) if you don't already have it
```R
install_github("ctlab/fgsea")
```

Enrichment analysis
```R
# download mouse symbol hallmarks gene sets
download.file("https://s3.amazonaws.com/cnobjects/singleCellNet/resources/mouse_symbols_H_v5p2_Dec_11_2017.rda")

gsHallmarks<-utils_loadObject("mouse_symbols_H_v5p2_Dec_11_2017.rda")

# limit to genes included in both the expression data and in the hallmarks lists
gsHall<-lapply(gsHallmarks, intersect, rownames(xdiff[[1]]))

system.time(enrPatt<-fgsea.wrapper.set(xdiff, gsHall, minSize=20, nPerm=1e3))
  user  system elapsed 
 34.546   0.814   6.781 

esMatPat<-ks.extract.more(enrPatt, sigThresh=1.1, sigType='padj', gsColName='pathway', esColName='NES')

# heatmap the enrichment scores, color only those that are significant
hm_enr(esMatPat, 0.05, cRows=T, cCols=T, fsr=6)
```
<img src="md_img/hm_fgsea_tm.png">



