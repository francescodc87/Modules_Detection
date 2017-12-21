
load("info_clusters.Rdata")
load("Modules.Rdata")
load("modules_adj_pvalues.Rdata")

N <- nrow(Modules)
Nclusters.in <- NULL
N.MiBiG.hits <- NULL
max.pvals <- NULL
ClassesofCompounds <- NULL
N.compound.classes <- NULL
ShannonsEntropy <- NULL

for(i in 1:N){
  tmp.module <- Modules[i,!is.na(Modules[i,])]
  ind <- 1:ncol(Clusters.complete)
  
  for(j in 1:length(tmp.module)){
    ind.tmp <- which(Clusters.complete==tmp.module[j], arr.ind = T)
    ind <- intersect(ind, ind.tmp[,2])
  }
  
  Nclusters.in <- c(Nclusters.in, length(ind))
  N.MiBiG.hits <- c(N.MiBiG.hits,
                    length(which(!is.na(Clusters.complete["MiBiG id",ind]))))  
  
  ind.pval <- which(rownames(coloc.pval)%in%tmp.module)
  pval.coloc <- coloc.pval[ind.pval, ind.pval]
  pval.coloc <- pval.coloc[upper.tri(pval.coloc)]
  pval.neigh <- neigh.pval[ind.pval, ind.pval]
  pval.neigh <- pval.neigh[upper.tri(pval.neigh)]
  pvals <- apply(cbind(pval.coloc, pval.neigh), 1, min, na.rm=T)
  max.pvals <- c(max.pvals, max(pvals, na.rm=T))
  Compound.classes <- Clusters.complete["Compound class", ind]
  Compound.classes <- unlist(strsplit(Compound.classes, split="-"))
  ClassesofCompounds <- c(ClassesofCompounds, paste(unique(Compound.classes), collapse = "-"))  
  N.compound.classes <- c(N.compound.classes, length(unique(Compound.classes)))
  tmp.freq <- table(Compound.classes)/length(Compound.classes)
  ShannonsEntropy <- c(ShannonsEntropy, sum(-tmp.freq*log(tmp.freq)))
    
}


save(len.mods, Nclusters.in,
           N.MiBiG.hits, max.pvals, 
           ClassesofCompounds,
           N.compound.classes,
           ShannonsEntropy,
           file="Metric_part.Rdata")



#### adding_composition_ranking
load("Metric_part.Rdata")

COG.ANNOTATIONs <- read.csv("COG_annotation_final.csv")
Classes.COGS <- levels(COG.ANNOTATIONs[,2])

adding <- matrix(0,nrow(Modules), length(Classes.COGS))
colnames(adding) <- Classes.COGS


for(i in 1:nrow(Modules)){
  tmp.module <- Modules[i,!is.na(Modules[i,])]
  ind.p <- which(COG.ANNOTATIONs[,1]%in%tmp.module)
  
  for(j in ind.p){
    ind.c <- which(colnames(adding)==COG.ANNOTATIONs[j,2])
    adding[i,ind.c]<- adding[i,ind.c]+1
  }
  adding[i,] <- (adding[i,]/len.mods[i])*100
  if(!i %% 1000){ cat("\n", (i/dim(adding)[1])*100, "%")}
}


save(adding, file="cog_annotation_metrics.Rdata")



Metrics.Modules <- cbind(len.mods, ShannonsEntropy, Nclusters.in,
                         N.MiBiG.hits,
                         max.pvals,
                         N.compound.classes,
                         ClassesofCompounds, adding)

colnames(Metrics.Modules) <- c("length", "Shannon's Entropy", "N cluster hit",
                               "N clusters MIBIG", "max pvalues", "N compound classes",
                               "Compound classes", colnames(adding))

save(Modules, Metrics.Modules, file="all_detected_modules.Rdata")

load("all_detected_modules.Rdata")

ind.keep <- which(Nclusters.in>1)

Modules <- Modules[ind.keep,]
Metrics.Modules <- Metrics.Modules[ind.keep,]

nrow(Modules)
Modules.names <- paste("Module",1:nrow(Modules), sep="_")

rownames(Modules) <- Modules.names
rownames(Metrics.Modules) <- Modules.names

save(Modules, Metrics.Modules, file="filtered_modules_2hits.Rdata")


###creating Ranks matrix for the MIB score









