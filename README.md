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
### trimming
this step starts from the file clustersCOGSnr.csv and creates the clustersCOGSnr_trimmed.csv where:
1- all clusters with less than 2 different cogs are removed
2- if 2 or more clusters have the same smCOG composition, only the shortest is kept
3- when the same smCOG is repeated subsequentially more than once in a cluster
    it is replaced with the empty cell "-" and concatenated at the end of the
    cluster as c(cluster, "-", "smCOG").
    Even if we are losing some neighbouring interactions and slightly
    modifying the cluster's topology, this makes the the p-value evaluation
    extremely simpler and extremely faster

###Example
 - Item 1
 - Item 2
  - Sub Item 1
  - Sub Item 2

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
