# Modules Detection
Here you can find all the R code and data used fo the module detection algorithm
described in "Computational identification of co-evolving multi-gene modules in microbial biosynthetic gene clusters".
Before starting, it is necessary to download everything contained here.

## Prerequisites
Several packages needs to be installed and loaded in order to be able to use the code here described.
* [parallel](https://cran.r-project.org/web/views/HighPerformanceComputing.html)
* [foreach](https://cran.r-project.org/web/packages/foreach/index.html)
* [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)
* [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html)
* [igraph](http://igraph.org/r/)

## Getting Started
Keep in mind that all intermediate files are saved, so it is
not necessary to replicate the whole analysis. 

First it is necessary to load all the necessary packages:

```
library(parallel)
library(foreach)
library(doParallel)
library(Matrix)
library(igraph)
```
After, one should navigate to the folder where all the data and the scripts are stored.
```
setwd("/path/to/chosen/folder/")
```
### Trimming
This step starts from the file containing all the BGCs annotated according smCOGs (clustersCOGSnr.csv) and creates the clustersCOGSnr_trimmed.csv where:
1. all clusters with less than 2 different cogs are removed;
1. if 2 or more clusters have the same smCOG composition, only the shortest is kept;
1. when the same smCOG is repeated subsequentially more than once in a cluster it is replaced with the empty cell "-" and concatenated at the end of the cluster as c(cluster, "-", "smCOG"). Even if we are losing some neighbouring interactions and slightly modifying the cluster's topology, this makes the the p-value evaluation extremely simpler and extremely faster.

To perform this step, one just need to paste the following:

```
source("trimming.R")
rm(list = ls())
```
### Counting colocalization interactions
This step takes as input the clustersCOGSnr_trimmed.csv file, counts all the coloc interactions between smCOGs and puts them in a matrix. The results are saved both in coloc_count.Rdata and coloc_count.csv.

```
source("counting_coloc.R")
rm(list = ls())
```
### Counting adjacency interactions
This step takes as input the clustersCOGSnr_trimmed.csv file, counts all the neigh interactions between smCOGs and puts them in a matrix. The results are saved both in neigh_count.Rdata and neigh_count.csv.
```
source("counting_neigh.R")
rm(list = ls())
```

### Computing colocalization p-values
This step takes as input the clustersCOGSnr_trimmed.csv and coloc_count.Rdata files and computes the colocalization p-values in both directions. The computed p-values are saved into coloc_pval.Rdata and coloc_pval.csv.
```
source("coloc_pval.R")
rm(list = ls())
```
# NOT FINISHED YET


## Built With
* [R](https://www.r-project.org/)
* [RStudio](https://www.rstudio.com/)


## Authors
* **Francesco Del Carratore**
* **Konrad Zych**
* **Matthew Cummings**
* **Eriko Takano**
* **Marnix Medema**
* **Rainer Breitling**
