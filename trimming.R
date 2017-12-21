### here I am trimming the cluster file in a different way.
clusters <- as.matrix(read.csv("clustersCOGSnr.csv", sep=";",
                               stringsAsFactors = F, header = F, row.names = 1, na.strings = "" ))



#### First I have to get rid of redundancy, the idea is to remove all the clusters
### that show the exact same composition in terms of defined clusters (ignoring -)
tmp <- clusters[1,]
ind <- which(!is.na(tmp) & tmp!= "-")
size <- length(which(!is.na(tmp)))
tmp <- tmp[ind]
tmp <- unique(tmp)
tmp <- sort(tmp)
Composition <-tmp


for(i in 2:dim(clusters)[1]){
    tmp <- clusters[i,]
    ind <- which(!is.na(tmp) & tmp!= "-")
    size <- c(size,length(which(!is.na(tmp))))
    tmp <- tmp[ind]
    tmp <- unique(tmp)
    tmp <- sort(tmp)
    Composition <- qpcR:::rbind.na(Composition, tmp) ### binding while adding NAs
}
rownames(Composition) <- rownames(clusters)
Comp <- unique(Composition, MARGIN = 1)
Composition[is.na(Composition)] <- ""
Comp[is.na(Comp)] <- ""

dim(Composition)
dim(Comp)
### I should end up with 12960 clusters (from 14598)


keep <- NULL
for(j in 1:dim(Comp)[1]){
     TMP <- which(apply(Composition, 1, function(x) all(x==Comp[j,], na.rm = T)))
     if(length(TMP)==1) {
         keep <- c(keep, TMP)
     }else{
         tmp <- which(size[TMP]==min(size[TMP]))
         keep <- c(keep, TMP[tmp[1]])
     }
 }
 keep <- sort(keep)
 
 
 dim(clusters)
 clusters.trimmed <- clusters[keep,]
 dim(clusters.trimmed)


###### When I arrive here I removed all the clusters that show the same composition (I kept the shorter one, i.e. the one with
###### less repetitions)
 
 
###  Now I remove all the clusters that show less than 2 COGs.
 
 kept <- NULL
 
 for (i in 1: dim(clusters.trimmed)[1]){
     tmp <- clusters.trimmed[i,]
     ind <- which(!is.na(tmp))
     tmp <- unique(tmp[ind])
     check <- length(which(tmp!="-"))
     if (check >= 2){
         kept <- c(kept, i)
     }
 }
 
 dim(clusters.trimmed)
 clusters.trimmed <- clusters.trimmed[kept,]
 dim(clusters.trimmed)

 
  
 clusters.trimmed.old <- clusters.trimmed


###### Now I have to check each cluster for repetiontions... 

Nclust <- dim(clusters.trimmed)[1]
clusters.trimmed <- cbind(clusters.trimmed, matrix(NA,Nclust,100))
Lclust <- dim(clusters.trimmed)[2]


for(j in 1:Nclust){
    clust <- clusters.trimmed[j,]
    ind <- which(!is.na(clust))
    clust <- clust[ind]
    ind.full <- which(clust!="-")
    clust2 <-clust[ind.full]
    repeated <- clust2[which(duplicated(clust2))]
    repeated <- unique(repeated)
    for(k in 1:length(repeated)){
        ind2 <- which(clust==repeated[k])
        if(length(ind2)==2 & ind2[2]==(ind2[1]+1)){
            clust[ind2] <- c("r",repeated[k])
            keep <- which(clust!="r")
            clust <- clust[keep]
        }else if(length(ind2)!=0){
            clust[ind2]<-"-"
            clust <- c(clust,c("-", repeated[k]))
        }
    }
    size <- length(clust)
    clust <- c(clust, rep(NA,Lclust-size))
    clusters.trimmed[j,] <- clust    
}

i <- 1
check <- FALSE
while(!check){
    check <- all(is.na(clusters.trimmed[,i]))
    i <- i+1
}
clusters.trimmed <- clusters.trimmed[,1:(i-2)]




## identifing the smCOGs left
ind <- which(!is.na(clusters.trimmed) & clusters.trimmed!= "-")
smCOGs <- sort(unique(as.numeric(clusters.trimmed[ind])))
smCOGsNames <- NULL

for (i in 1:length(smCOGs)){
    tmp <-as.character(10000+smCOGs[i])
    smCOGsNames <- c(smCOGsNames, paste("SMCOG", tmp, sep = ""))
}

names(smCOGs)<- smCOGsNames
rm(smCOGsNames, ind, i)


###Now I can save
write.csv(clusters.trimmed, "clustersCOGSnr_trimmed.csv")
rm(list = ls())


