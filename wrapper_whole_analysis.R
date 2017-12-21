#### this file calls all the scripts in order to replicate the module analysis.
#### Keep in mind that I am also saving all intermediate files, so it is not
#### not necessary to replicate the whole analysis 

# list of packages needed
#library(parallel)
#library(foreach)
#library(doParallel)
#library(Matrix)
#library(igraph)




### here we need to set up the folder where the needed files are stored
### and the results saved

setwd("/home/mqbpqfd2/Desktop/Modules_scripts_arranged/sharing_Modules_analysis/")


### trimming ###
# this code starts from the file clustersCOGSnr.csv and creates the
# clustersCOGSnr_trimmed.csv where:
# 1- all clusters with less than 2 different cogs are removed
# 2- if 2 or more clusters have the same smCOG composition, only the shortest is kept
# 3- when the same smCOG is repeated subsequentially more than once in a cluster
#    it is replaced with the empty cell "-" and concatenated at the end of the
#    cluster as c(cluster, "-", "smCOG").
#    Even if we are losing some neighbouring interactions and slightly
#    modifying the cluster's topology, this makes the the p-value evaluation
#    extremely simpler and extremely faster

source("trimming.R")
rm(list = ls())


####  counting coloc interactions
## this script takes as input the clustersCOGSnr_trimmed.csv file,
### counts all the coloc interactions between smCOGs and puts them in a matrix
### the results are saved both in coloc_count.Rdata and coloc_count.csv
source("counting_coloc.R")
rm(list = ls())



####  counting neigh interactions
## this script takes as input the clustersCOGSnr_trimmed.csv file,
### counts all the neigh interactions between smCOGs and puts them in a matrix
### the results are saved both in neigh_count.Rdata and neigh_count.csv
source("counting_neigh.R")
rm(list = ls())


#### computing coloc p-values
## this script takes as input the clustersCOGSnr_trimmed.csv and
## coloc_count.Rdata files and computes the coloc p-values in both directions
## and saves them into coloc_pval.Rdata and coloc_pval.csv
source("coloc_pval.R")
rm(list = ls())



#### computing neigh p-values
## this script takes as input the clustersCOGSnr_trimmed.csv and
## neigh_count.Rdata files and computes the coloc p-values in both directions
## and saves them into neigh_pval.Rdata and neigh_pval.csv
source("neigh_pval.R")
rm(list = ls())


#####handling coloc p-values
# this script takes as input the coloc_pval.Rdata file. First, it creates a simmetric matix
# of p-values where the (i,j) element is subtituted by the max((i,j),(j,i)) and saves it into
# coloc_pval_MAX.Rdata and coloc_pval_MAX.csv. Second, it takes the last matrix and applies the
# BY correction to the p-values. This is saved in coloc_pval_adj.Rdata and coloc_pval_adj.csv
source("handling_coloc_pvalues.R")
rm(list = ls())


#####handling neigh p-values
# this script takes as input the neigh_pval.Rdata file. First, it creates a simmetric matix
# of p-values where the (i,j) element is subtituted by the max((i,j),(j,i)) and saves it into
# neigh_pval_MAX.Rdata and neigh_pval_MAX.csv. Second, it takes the last matrix and applies the
# BY correction to the p-values. This is saved in neigh_pval_adj.Rdata and neigh_pval_adj.csv
source("handling_coloc_pvalues.R")
rm(list = ls())


###putting the adj p-values matrices in the same file

load("coloc_pval_adj.Rdata")
load("neigh_pval_adj.Rdata")
save(coloc.pval,neigh.pval, file = "modules_adj_pvalues.Rdata")
rm(list = ls())


#### detecting modules
# from modules_adj_pvalues.Rdata to Modules.Rdata
## here the modules detection happens. Starting from an arbitrary
## threshold (.1), a simmetric binary matrix is built where the i,j element
## is 1 if at least one of the two pvalues (neigh or coloc) is lower than threshold.
## This binary matrix is than used to create a graph. The the igraph package is
## used to find all the maximal cliques (fully connected sub-graphs there are not
## subsets of other fully connected sub-graphs). These cliques are modules!!!
## this procedure is than repeated with every possible threshold < .1
## At the end 197564 putative modules are found.

source('detect_modules_function.R')
source("detect_modules.R")
rm(list = ls())


#### computing Metrics
# this scripts takes as inputs: info_clusters.Rdata, Modules.Rdata,
# modules_adj_pvalues.Rdata, COG_annotation_final.csv
# it computes a number of different metrics used to rank and filter the modules:
## length, Shannon's entropy, N cluster hit, N clusers MIBIG, max pval, N compound classes,
## % CORE cogs, % CORE/TAILORING cogs, % MIXED cogs, % OTHER cogs, % REGULATOR cogs,
## % REGULATOR/TAILORING cogs, % TAILORING cogs, % TAILORING/CORE cogs, % TAILORING/REGULATOR cogs,
## % TAILORING/TRANSPORT cogs and % TRANSPORT cogs
##The file all_detected_modules.Rdata contains both module composition and metrics per each
##of the 197564 detected modules. The filtered_modules_2hits.Rdata instead contains all modules
## found in at least 2 clusters (they are modules otherwise)
source("Computing_Metrics.R")
rm(list = ls())
### check if the metrics are considering the very last MiBiG mapping...now it's wrong


#### computing MIB score
# this score it basically is a weighted rank sum
# first I  have to compute the ranks for 
source("compute_ranks.R")
rm(list = ls())

#define the weights
       
weights <- c(2, #"length"
             15, #"Shannon's Entropy"
             10, #"N cluster hit"
             5, #"max pvalues"
             0, #"N compound classes"
             0, # "CORE"
             0, # "CORE/TAILORING"
             0, #"MIXED"
             0,  #"OTHER"
             0, #"REGULATOR"
             0, # "REGULATOR/TAILORING"
             10, #"TAILORING"
             0, #"TAILORING/CORE"
             0, #"TAILORING/REGULATOR"
             0, #"TAILORING/TRANSPORT"
             0) #"TRANSPORT"   

#actual MIB score - final output
source("MIB_score.R")




load("final_results_modules.Rdata")




BP <- cbind(as.numeric(Metrics.Modules[,2]),
            as.numeric(Metrics.Modules[,3]),
            as.numeric(Metrics.Modules[,4]),
            as.numeric(Metrics.Modules[,5]),
            as.numeric(Metrics.Modules[,6]),
            as.numeric(Metrics.Modules[,7]))
boxplot(BP)




