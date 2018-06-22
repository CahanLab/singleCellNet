# singleCellNet

### Introduction
See [CellNet](https://github.com/pcahan1/CellNet) for an introduction to CellNet, how to use it on bulk RNA-Seq, and how to analyze single cell RNA-Seq (scRNA-Seq) data with classifiers trained on bulk RNA-Seq. Here, we illustrate

- how to build and assess single cell classifiers

- how to use these classifiers to quantify 'cell identity' from query scRNA-Seq data

- how to cluster scRNA-Seq data using our 'cluster by competition' method

### DATA

In this example, we use a subset of the Tabula Muris data to train singleCellNet. To learn more about the Tabula Muris project, see the [manuscript])(https://www.biorxiv.org/content/early/2018/03/29/237446). As query data, we use scRNA-Seq of kidney cells as reported in [Park et al 2018](https://www.ncbi.nlm.nih.gov/pubmed/29622724). You can download this data here:

| APPLICATION | METADATA | EXPRESSION |
|-------------------------------------|
| Query       | [metadata](https://s3.amazonaws.com/cnobjects/singleCellNet/examples/sampTab_Park_MouseKidney_062118.rda) | [expression data](https://s3.amazonaws.com/cnobjects/singleCellNet/examples/GSE107585_Mouse_kidney_single_cell_datamatrix.txt.gz) |
| Training    | [metadata](https://s3.amazonaws.com/cnobjects/singleCellNet/examples/sampTab_TM_053018.rda) | [expression data](https://s3.amazonaws.com/cnobjects/singleCellNet/examples/expTM_Raw_053018.rda) |

N.B. The query expression data needs to be decompressed before loading it into R. 


#### Setup
```R
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





