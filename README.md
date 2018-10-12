# singleCellNet

### Introduction
See [CellNet](https://github.com/pcahan1/CellNet) for an introduction to CellNet, how to use it on bulk RNA-Seq, and how to analyze single cell RNA-Seq (scRNA-Seq) data with classifiers trained on bulk RNA-Seq. Here, we illustrate

- how to build and assess single cell classifiers 

- how to build and assess cross-species single cell classifiers

- how to use these classifiers to quantify 'cell identity' from query scRNA-Seq data

- how to cluster scRNA-Seq data using our 'cluster by competition' method

### DATA

In this example, we use a subset of the Tabula Muris data to train singleCellNet. To learn more about the Tabula Muris project, see the [manuscript])(https://www.biorxiv.org/content/early/2018/03/29/237446). As query data, we use scRNA-Seq of kidney cells as reported in [Park et al 2018](https://www.ncbi.nlm.nih.gov/pubmed/29622724). We also provide an example of classifying human, bead enriched PBMCs (from https://www.ncbi.nlm.nih.gov/pubmed/28091601). You can download this data here:

| APPLICATION | METADATA | EXPRESSION |
|-------------|----------|------------|
| Query       | [metadata](https://s3.amazonaws.com/cnobjects/singleCellNet/examples/sampTab_Park_MouseKidney_062118.rda) | [expression data](https://s3.amazonaws.com/cnobjects/singleCellNet/examples/expDat_Park_MouseKidney_062218.rda") |
| Training    | [metadata](https://s3.amazonaws.com/cnobjects/singleCellNet/examples/sampTab_TM_053018.rda) | [expression data](https://s3.amazonaws.com/cnobjects/singleCellNet/examples/expTM_Raw_053018.rda) |
| cross-species | [human-mouse orthologs](https://s3.amazonaws.com/cnobjects/singleCellNet/examples/human_mouse_genes_Jul_24_2018.rda)| Query (human bead-purified PBMC from 10x) | [metadata](https://s3.amazonaws.com/cnobjects/singleCellNet/examples/stDat_beads_mar22.rda) | [expression data](https://s3.amazonaws.com/cnobjects/singleCellNet/examples/6k_beadpurfied_raw.rda) |

#### Setup
```R
library(fgsea)
library(devtools)
install_github("pcahan1/singleCellNet", ref="master", auth="your_token")
library(singleCellNet)

library(RColorBrewer)
library(pheatmap)
library(randomForest)
library(viridis)
library(ggplot2)
library(dplyr)
library(pROC)
library(viridis)
library(patchwork)
library(DescTools)
library(Matrix)
library(parallel)

mydate<-utils_myDate()
```

#### Fetch the data if you have not already done so
```R
download.file("https://s3.amazonaws.com/cnobjects/singleCellNet/examples/sampTab_Park_MouseKidney_062118.rda", "sampTab_Park_MouseKidney_062118.rda")

download.file("https://s3.amazonaws.com/cnobjects/singleCellNet/examples/expMatrix_Park_MouseKidney_Oct_12_2018.rda", "expMatrix_Park_MouseKidney_Oct_12_2018.rda")

download.file("https://s3.amazonaws.com/cnobjects/singleCellNet/examples/expMatrix_TM_Raw_Oct_12_2018.rda", "expMatrix_TM_Raw_Oct_12_2018.rda")

download.file("https://s3.amazonaws.com/cnobjects/singleCellNet/examples/sampTab_TM_053018.rda", "sampTab_TM_053018.rda")

## For cross-species analyis:
download.file("https://s3.amazonaws.com/cnobjects/singleCellNet/examples/human_mouse_genes_Jul_24_2018.rda", "human_mouse_genes_Jul_24_2018.rda")

download.file("https://s3.amazonaws.com/cnobjects/singleCellNet/examples/6k_beadpurfied_raw.rda", "6k_beadpurfied_raw.rda")

download.file("https://s3.amazonaws.com/cnobjects/singleCellNet/examples/stDat_beads_mar22.rda", "stDat_beads_mar22.rda")

```

#### Load query data
```R
stPark<-utils_loadObject("sampTab_Park_MouseKidney_062118.rda")
expPark<-utils_loadObject("expMatrix_Park_MouseKidney_Oct_12_2018.rda")
dim(expPark)
[1] 16272 43745

genesPark<-rownames(expPark)
```

#### Load the training data
```R
expTMraw<-utils_loadObject("expMatrix_TM_Raw_Oct_12_2018.rda")
dim(expTMraw)
[1] 23433 24936

stTM<-utils_loadObject("sampTab_TM_053018.rda")
dim(stTM)
[1] 24936    17

stTM<-droplevels(stTM)
```

#### Find genes in common to the data sets and limit analysis to these
```R
commonGenes<-intersect(rownames(expTMraw), genesPark)
length(commonGenes)
[1] 13831


expPark<-expPark[commonGenes,]
expTMraw<-expTMraw[commonGenes,]
```

#### Split for training and assessment, and transform training data
```R
stList<-splitCommon(stTM, ncells=100, dLevel="newAnn")
stTrain<-stList[[1]]
expTrain<-expTMraw[,rownames(stTrain)]


system.time(tmpX<-weighted_down(expTrain, 1.5e3, dThresh=0.25))
   user  system elapsed 
  4.837   0.845   5.711

system.time(expTrain<-trans_prop(tmpX, 1e4))
   user  system elapsed 
  1.486   0.645   2.136
```

#### Find the best set of classifier genes
```R

system.time(cgenes2<-findClassyGenes(expTrain, stTrain, "newAnn", topX=10))
  user  system elapsed 
 38.721  10.045  48.767

cgenesA<-cgenes2[['cgenes']]
grps<-cgenes2[['grps']]
length(cgenesA)
[1] 476

# heatmap these genes
hm_gpa_sel(expTrain, cgenesA, grps, maxPerGrp=5, toScale=T, cRow=F, cCol=F,font=4)
```
<img src="md_img/hm_tabulaMuris.png">

#### Find the best pairs
```R
expT<-as.matrix(expTrain[cgenesA,])
dim(expT)
[1]  476 3036

system.time(xpairs<-ptGetTop(expT, grps, topX=25, sliceSize=5000))
    user   system  elapsed 
1671.187 1406.671  154.199

length(xpairs)
[1] 799
```

#### TSP transform the training data
```R
system.time(pdTrain<-query_transform(expT[cgenesA, ], xpairs))

dim(pdTrain)
[1]  799 3036

 ```

#### Train the classifier
```R
system.time(rf_tspAll<-sc_makeClassifier(pdTrain[xpairs,], genes=xpairs, groups=grps, nRand=100, ntrees=1000))
  user  system elapsed 
166.643   0.248 166.866
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


#### Attribution plot
```R
plot_attr(classRes_val_all, stTest, nrand=nrand, dLevel="newAnn", sid="cell")
```
<img src="md_img/attribution_val_101218.png">

#### UMAP by category
```R
system.time(umPrep<-prep_umap_class(classRes_val_all, stTest, nrand=nrand, dLevel="newAnn", sid="cell", topPC=5))
user  system elapsed 
109.500   3.588 113.067 

plot_umap(umPrep)
```
<img src="md_img/umap_val_101218.png">


#### Assess classifier
```R
newSampTab<-makeSampleTable(classRes_val_all, stTest, nrand, "cell")
tm_heldoutassessment <- assessmentReport_comm(classRes_val_all, newSampTab, classLevels='newAnn',dLevelSID='cell')
plot_PRs(tm_heldoutassessment)
```
<img src="md_img/pr_101218.png">

```R
plot_metrics(tm_heldoutassessment, method = "tsp_rf", ylimForMultiLogLoss = 1000)
```
<img src="md_img/metrics_101218.png">


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

### Cross-species classification

Load the human query data
```R
stQuery<-utils_loadObject("stDat_beads_mar22.rda")
expQuery<-utils_loadObject("6k_beadpurfied_raw.rda") # use Matrix if RAM low
dim(expQuery)
[1] 32643  6000

expTMraw<-utils_loadObject("expMatrix_TM_Raw_Oct_12_2018.rda") # reload training

```

Load the ortholog table and convert human gene names to mouse ortholog names, and limit analysis to genes in common between the training and query data.
```R
oTab<-utils_loadObject("human_mouse_genes_Jul_24_2018.rda")
dim(oTab)
[1] 16688     3

aa = csRenameOrth(expQuery, expRawTM, oTab)
expQuery <- aa[['expQuery']]
expTrain <- aa[['expTrain']]
```

Limit anlaysis to a subset of the TM cell types
```R
cts<-c("B cell",  "cardiac muscle cell", "endothelial cell", "erythroblast", "granulocyte", "hematopoietic precursor cell", "late pro-B cell", "limb_mesenchymal", "macrophage", "mammary_basal_cell", "monocyte", "natural killer cell", "T cell", "trachea_epithelial", "trachea_mesenchymal")

stTM2<-filter(stTM, newAnn %in% cts)
stTM2<-droplevels(stTM2)
rownames(stTM2)<-as.vector(stTM2$cell) # filter strips rownames

expTrain<-expTrain[,rownames(stTM2)]
dim(expTrain)
[1] 14550 15161
```

Split into training and validation, normalize training data, and find classy genes
```R
stList<-splitCommon(stTM2, ncells=100, dLevel="newAnn")
stTrain<-stList[[1]]
dim(stTrain)

[1] 1457   17

expTMnorm<-trans_prop(weighted_down(expTrain[,rownames(stTrain)], 1.5e3, dThresh=0.25), 1e4)


system.time(cgenes2<-findClassyGenes(expTMnorm, stTrain, "newAnn", topX=10))
   user  system elapsed 
 12.735   3.114  15.847 


cgenesA<-cgenes2[['cgenes']]
grps<-cgenes2[['grps']]
length(cgenesA)
[1] 244

```


find best pairs and transform query data, and train classifier
```R
system.time(xpairs<-ptGetTop(expTMnorm[cgenesA,], grps, topX=25, sliceSize=5000))
   user  system elapsed 
117.248 120.292 115.127 

pdTrain<-query_transform(expTrain[cgenesA, rownames(stTrain)], xpairs)

dim(pdTrain)
[1]  375 1457

nrand = 50
system.time(rf_tspAll<-sc_makeClassifier(pdTrain[xpairs,], genes=xpairs, groups=grps, nRand=50, ntrees=1000))
   user  system elapsed 
 18.321   0.057  18.373
 ```

Apply to held out data
```R
stTest<-stList[[2]]

system.time(expQtransAll<-query_transform(expTrain[cgenesA,rownames(stTest)], xpairs))
   user  system elapsed 
  2.744   0.061   2.806 

nrand<-50
system.time(classRes_val_all<-rf_classPredict(rf_tspAll, expQtransAll, numRand=nrand))
   user  system elapsed 
  8.015   0.178   8.191 


sla<-as.vector(stTest$newAnn)
names(sla)<-rownames(stTest)
slaRand<-rep("rand", nrand)
names(slaRand)<-paste("rand_", 1:nrand, sep='')

sla<-append(sla, slaRand)

# heatmap classification result
sc_hmClass(classRes_val_all, sla, max=300, font=7, isBig=TRUE)
```
<img src="md_img/hmClass_CS_heldOut_101218.png">

Apply to human query data
```R
system.time(expQueryTrans<-query_transform(expQuery[cgenesA,], xpairs))
  user  system elapsed 
  0.149   0.027   0.176 
  
nqRand<-50
system.time(crHS<-rf_classPredict(rf_tspAll, expQueryTrans, numRand=nqRand))
   user  system elapsed 
  2.390   0.068   2.456 

# heatmap classification result
sgrp<-as.vector(stQuery$prefix)
names(sgrp)<-rownames(stQuery)
grpRand<-rep("rand", nqRand)
names(grpRand)<-paste("rand_", 1:nqRand, sep='')
sgrp<-append(sgrp, grpRand)
sc_hmClass(crHS, sgrp, max=5000, isBig=TRUE, cCol=F, font=8)
```
<img src="md_img/hmClass_CS_090618.png">

Note that the macrophage category seems to be promiscuous in the mouse held out data, too.



