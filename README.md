# singleCellNet

### Introduction
See [CellNet](https://github.com/pcahan1/CellNet) for an introduction and bulk RNA-Seq version. Here, we illustrate how to make single cell classifiers.

#### Processing pipeline

1. Load
2. Wash 
  * down sample
  * exclude undetected genes
  * apply selected transform
    * proportional (default)
    * rank
    * zscore
    * binarize
  * results in object that has
    * arg list
    * expMatrix
    * transMethod
3. Chop (dimension reduction)
  * methods
    * PCA
    * tSNE
  * results in a list of
    * choppedData
    * arg lis
    * varGenes
4. Steam (assign cells to groups)
  * methods
    * manual
    * mclust
    * dynamic tree cut
  * results in a list of
    * sample table
    * arg list
    * optimized parameters (mclust=G and model shape, dt=minModSize and deepSplit, dbscan=eps and minPts)
    * method name
5. Butter (make and assess classifiers)
6. Toss (classify new samples)
7. Mix (integrate new training data)


#### Setup

```R
    install_github("pcahan1/singleCellNet", ref="master")
    library(singleCellNet)
    library(cellrangerRkit)
    library(Rtsne)
    library(ggplot2)
    library(pheatmap)
    library(dbscan)
    library(RColorBrewer)
    library(WGCNA)
    library(mclust)
    library(randomForest)
```
    
#### Load data
```R
    rawDat<-mergeLoad10x("pathTo/10x_public/Zheng_2016/bead_purified/", c("bcell_cd19", "cd34", "monocytes_cd14", "nkcell_cd56", "tcell_cd4_helper", "tcell_cd8_cytotoxic"), nCells=1e3))
```
#### Basic transform: normalize to total counts, then Log(1+norm), after down sampling
```R
    expDat<-rawDat[['expDat']]
    stDat<-rawDat[['sampTab']]

    expDn<-weighted_down(expDat,1e3)

    expNorm<-trans_dNorm(expDn)
    expNorm<-log10(1+expNorm)
```

#### Wash
```R
    pwashed<-prewash(expNorm, stDat)
    washedDat<-wash(expNorm, pwashed)
```

#### Chop, Steam, and assess classifiers based only on expProp
```R
    cAssAll<-pipe_cAss(washedDat, stDat)
    ggplot(cAssAll, aes(x=group, y=classDiff)) + geom_boxplot(alpha=.75,colour="black", outlier.colour="#BDBDBD", outlier.size=.5) + xlab("cluster") + theme_bw() + facet_wrap( ~ method)
```

#### Chop and Steam, useful when assessing various wash methods
```R
    steamed<-pipe_steam_list(washedDat, stDat, topPC=20)
```

#### make classifiers and assess -- expProp
```R
    classAssProp<-pipe_cAss_all(steamed, expDat, stDat)
    ggplot(classAssProp, aes(x=method, y=classDiff)) + geom_boxplot(alpha=.75,colour="black", outlier.colour="#BDBDBD", outlier.size=.5) + xlab("cluster") + theme_bw()
```
#### Binary data
```R
    expBin<-trans_binarize(expDn, threshold=1)
    classAssBin<-pipe_cAss_all(steamed, expBin, stDat)
```

#### zscore data
```R
    expZscore<-trans_zscore(expDn)
    # make sure gene names are here -- need to fix this
    rownames(expZscore)<-rownames(expDn)
    classAssZscore<-pipe_cAss_all(steamed, expZscore, stDat)
````

#### compare all methods
```R
    cAssBound<-cbind(classAssProp, wash=rep("prop", nrow(classAssProp)))
    cAssBound<-rbind(cAssBound, cbind(classAssBin, wash=rep("binary", nrow(classAssBin))))
    cAssBound<-rbind(cAssBound, cbind(classAssZscore, wash=rep("zscore", nrow(classAssZscore))))

    ggplot(cAssBound, aes(x=method, y=classDiff,fill=wash )) + geom_boxplot(alpha=.75,colour="black", outlier.colour="#BDBDBD", outlier.size=.5) + xlab("Steam method") + theme_bw()
```

Added on 06-01-17

Instructions for loading package, making a dbscan/binary data classifier, and applying it to a new data set

![](md_img/hm1.jpg)

![](md_img/hm2.jpg)


