---
title: "DMCHMM: Differentially Methylated CpG using Hidden Markov Model"
author: "Farhad Shokoohi"
date: "2020-03-07"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteEngine{knitr::knitr}
  %\VignetteIndexEntry{Sending Messages With Gmailr}
  %\usepackage[utf8]{inputenc}
---

DNA methylation studies have increased in number over the past decade thanks 
to the recent advances in next-generation sequencing (NGS) and microarray 
technology (MA), providing many data sets at high resolution, enabling 
researchers to understand methylation patterns and their regulatory roles in 
biological processes and diseases.  
Notwithstanding that diverse methods and software have created ample 
opportunities for researchers to do quantitative analysis, they make it 
difficult for practitioners to choose the one that is suitable and efficient 
in analyzing DNA methylation data. 
Having examined most of differentially methylation identification tools for 
bisulfite sequencing (BS-Seq) data, we observed several drawbacks in the 
existing analytic tools. To address these issues we have developed a novel 
differentially methylated CpG site identification tool which is based on Hidden
Markov models (HMM) called `DMCHMM`. This vignette provides some guidelines on 
how to use the package and analyze BS-Seq data. 

Following topics will be discussed in this vignette:

- **Reading BS-Seq data**
- **Simulating BS-Seq data**
- **Predicting methylation levels using HMM and EM algorithm**
- **Predicting methylation levels using HMM and MCMC algorithm**
- **Identifying DMCs**

## S4 Classes

Two different classes are defined by extending the 
`SummarizedExperiment-class`. The `BSData-class` is designed to hold BS-Seq 
data. Similarly, `cBSData-method` is defined to create a `BSData` object. 
This class includes two slots:
the `methReads`, a matrix with columns representing samples and rows 
representing genomic positions (CpG sites) and elements of matrix representing 
methylation counts at each position in each sample; 
the `totalReads`, a matrix with similar columns and rows except the elements 
representing total number of reads. 

## Reading BS-Seq data

For reading raw BS-Seq data we adopted The `readBismark` function from `BiSeq`
package. The `readBismark-method` reads samples stored in different files with 
six columns of *chromosome*, *start position*, *end position*, 
*methylation percentage*, *number of Cs* and *number of Ts*. 

Three data files are included in the `DMCHMM` package for illustration. 
The data can be imported using following code. 


```r
library(DMCHMM)
fn <- list.files(system.file("extdata",package = "DMCHMM"))
fn.f <- list.files(system.file("extdata",package="DMCHMM"), full.names=TRUE)
OBJ <- readBismark(fn.f, fn)
```

```
## Processing sample blk.BCU1568_BC_BS_1 ... 
## Processing sample blk.BCU173_TC_BS_1 ... 
## Processing sample blk.BCU551_Mono_BS_1 ... 
## Building BSData object.
```

```r
cdOBJ <- DataFrame(Cell = factor(c("BC", "TC","Mono"),
labels = c("BC", "TC", "Mono")), row.names = c("BCU1568","BCU173","BCU551"))
colData(OBJ) <- cdOBJ
OBJ
```

```
## class: BSData 
## dim: 25668 3 
## metadata(0):
## assays(2): totalReads methReads
## rownames(25668): 1 2 ... 25667 25668
## rowData names(0):
## colnames(3): BCU1568 BCU173 BCU551
## colData names(1): Cell
```

## Simulating BS-Seq data

The above data set only include one sample for each cell type. We need more 
samples to be able to compare their methylations and find DMCs. For illustration
we generate a sample of BS-Seq data as follows. 


```r
nr <- 150; nc <- 8
metht <- matrix(as.integer(runif(nr * nc, 0, 20)), nr)
methc <- matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
r1 <- GRanges(rep("chr1", nr), IRanges(1:nr, width=1), strand="*")
names(r1) <- 1:nr
cd1 <- DataFrame(Group=rep(c("G1","G2"),each=nc/2),row.names=LETTERS[1:nc])
OBJ1 <- cBSData(rowRanges=r1,methReads=methc,totalReads=metht,colData=cd1)
OBJ1
```

```
## class: BSData 
## dim: 150 8 
## metadata(0):
## assays(2): totalReads methReads
## rownames(150): 1 2 ... 149 150
## rowData names(0):
## colnames(8): A B ... G H
## colData names(1): Group
```


## Predicting methylation levels using HMM and EM algorithm

There are two approaches to smoothed the data before testing for DMCs. Either EM
or MCMC can be used to predict methylation levels utilizing HMM. The 
`methHMEM-method` which is developed to predict methylation levels. The output 
is a `BSDMCs-class` that can be either used to find DMCs or use MCMC algorithm 
to re-smooth the raw data. The process is as follows. 


```r
OBJ2 <- methHMEM(OBJ1, MaxK=2)
```

```
##   |                                                                              |                                                                      |   0%  |                                                                              |=========                                                             |  12%  |                                                                              |==================                                                    |  25%  |                                                                              |==========================                                            |  38%  |                                                                              |===================================                                   |  50%  |                                                                              |============================================                          |  62%  |                                                                              |====================================================                  |  75%  |                                                                              |=============================================================         |  88%  |                                                                              |======================================================================| 100%
```

```r
OBJ2
```

```
## class: BSDMCs 
## dim: 150 8 
## metadata(3): K Beta Pm
## assays(5): methReads totalReads methLevels methStates methVars
## rownames(150): 1 2 ... 149 150
## rowData names(0):
## colnames(8): A B ... G H
## colData names(1): Group
```

## Predicting methylation levels using HMM and MCMC algorithm

Although EM algorithm is a fast way to smooth the data but the results are not
as good as the MCMC algorithm. The MCMC algorithm, however, is slow. In order to
increase the speed, we first use `methHMEM-method` to find the HMM order for 
each sample and then we use `methHMCMC-method` to predict methylation levels. 
The procedure is as follows. 


```r
OBJ3 <- methHMMCMC(OBJ2)
```

```
##   |                                                                              |                                                                      |   0%  |                                                                              |=========                                                             |  12%  |                                                                              |==================                                                    |  25%  |                                                                              |==========================                                            |  38%  |                                                                              |===================================                                   |  50%  |                                                                              |============================================                          |  62%  |                                                                              |====================================================                  |  75%  |                                                                              |=============================================================         |  88%  |                                                                              |======================================================================| 100%
```

```r
OBJ3
```

```
## class: BSDMCs 
## dim: 150 8 
## metadata(3): K Beta Pm
## assays(5): methReads totalReads methLevels methStates methVars
## rownames(150): 1 2 ... 149 150
## rowData names(0):
## colnames(8): A B ... G H
## colData names(1): Group
```

## Identifying DMCs

Having smoothed the data using HMM, we run linear between predicted methylation
levels and grouping covariate. In case other covariates exist, one can use the 
`formula` argument to specify a linear model. When there is no covariates no 
action is required. The following command identifys the DMCs. The results are
stored in a `BSDMCs-class` and can be retrived by calling `metadata` command. 


```r
OBJ4 <- findDMCs(OBJ3)
```

```
##   |                                                                              |                                                                      |   0%  |                                                                              |=====                                                                 |   7%  |                                                                              |==========                                                            |  14%  |                                                                              |===============                                                       |  21%  |                                                                              |====================                                                  |  29%  |                                                                              |=========================                                             |  36%  |                                                                              |==============================                                        |  43%  |                                                                              |===================================                                   |  50%  |                                                                              |========================================                              |  57%  |                                                                              |=============================================                         |  64%  |                                                                              |==================================================                    |  71%  |                                                                              |=======================================================               |  79%  |                                                                              |============================================================          |  86%  |                                                                              |=================================================================     |  93%  |                                                                              |======================================================================| 100%
```

```
## Warning in fdrtool(x, statistic = "pvalue", plot = FALSE, verbose = FALSE):
## There may be too few input test statistics for reliable FDR calculations!
```

```r
head(metadata(OBJ4)$DMCHMM)
```

```
##   DMCs    pvalues   qvalues DMCsGroupG1vsG2 methDirGroupG1vsG2
## 1    0 0.80631622 0.8587848               0               hypo
## 2    0 0.38112884 0.7628627               0               hypo
## 3    0 0.02303307 0.4932287               0              hyper
## 4    0 0.12534059 0.7117866               0               hypo
## 5    0 0.43435342 0.7661639               0               hypo
## 6    0 0.13618504 0.7174863               0              hyper
```
