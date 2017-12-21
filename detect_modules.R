
library(Matrix)
library(igraph)

## alternative method to detect modules

load("modules_adj_pvalues.Rdata")

initial.threshold <- .1 ### this is chosen by me

Modules <- detect.modules(coloc.pval, neigh.pval, thr=initial.threshold)
### here I have 14094 
### I can do I complete look for sub-modules bigger than 3

len.mods <- NULL
for (i in 1:nrow(Modules)){
  len.mods <- c(len.mods, length(which(!is.na(Modules[i,]))))
}

save(Modules, len.mods, file="Modules_newApproach_initial.Rdata")

ind.to.explore  <- which(len.mods>3)
#### now for each of these I have to look for all the sub-modules I can find considering all
#### possible thresholds!

Modules.sub <- NULL
Modules.sub2 <-NULL
save(Modules.sub2, file="ModuleSub.Rdata")
flag<- 0

for(k in ind.to.explore){
  big.module <- Modules[k,!is.na(Modules[k,])]
  cat("\n", big.module, "\n", 
      (which(ind.to.explore==k)/length(ind.to.explore))*100, "%")
  ind.big <- which(rownames(coloc.pval)%in%big.module)
  
  coloc.pval.big <- coloc.pval[ind.big, ind.big]
  neigh.pval.big <- neigh.pval[ind.big, ind.big]
  
  thresholds <- unique(c(as.vector(coloc.pval.big[coloc.pval.big<.1]),
                               as.vector(neigh.pval.big[neigh.pval.big<.1])))
  thresholds <- thresholds[thresholds!=0]
  thresholds <- sort(thresholds, decreasing = T)
  
  
  for(i in 1:length(thresholds)){
  #  cat("\ncomputing with threshold", thresholds.known[i], i*100/length(thresholds.known), "%")
    tmp <- detect.modules(coloc.pval.big, neigh.pval.big, thresholds[i])
    #Modules.sub <- unique(rbind(Modules.sub,tmp))
    Modules.sub <- rbind(Modules.sub, tmp)
  }
  
  flag <- flag+1
  
  if(flag==20){
    load("ModuleSub.Rdata")
    Modules.sub2 <- rbind(Modules.sub, Modules.sub2)
    save(Modules.sub2, file="ModuleSub.Rdata")
    rm(Modules.sub2)
    Modules.sub <- NULL
    flag <- 0
  }
}

load("ModuleSub.Rdata")
Modules.sub2 <- rbind(Modules.sub, Modules.sub2)
save(Modules.sub2, file="ModuleSub.Rdata")
rm(Modules.sub2)

rm(tmp)


load("ModuleSub.Rdata")
Modules <- unique(rbind(Modules,Modules.sub2))
save(Modules, file="Modules.Rdata")
rm(Modules.sub2)
dim(Modules)
###I have to put togheter all the files here...after doing unique and stuff...
### next I have to prepare 


len.mods <- NULL
for(i in 1:nrow(Modules)){
  len.mods <- c(len.mods, length(which(!is.na(Modules[i,]))))
}

hist(len.mods)

save(Modules, len.mods,
     file="Modules.Rdata")




#### I could easily parallelise this for loop and make it extremely faster (for the future)








