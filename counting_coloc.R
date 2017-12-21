#defining the function
count.coloc.symmetric <- function(i,j,smCOGs,clusters){
    if(i==j){
        out <- 0
    }else{
        cogA <- smCOGs[i]
        cogB <- smCOGs[j]
        indA <- which(clusters==cogA, arr.ind = T)
        indB <- which(clusters==cogB, arr.ind = T)
        out <- length(intersect(indA[,1], indB[,1]))  
    }
    out
}



library(parallel)
library(foreach)
library(doParallel)
library(Matrix)
##### Counting colocalization interactions (asymmetric new)

clusters <- as.matrix(read.csv("clustersCOGSnr_trimmed.csv",
                               stringsAsFactors = F, header = T, row.names = 1 ))

ind <- which(!is.na(clusters) & clusters!= "-")
smCOGs <- sort(unique(as.numeric(clusters[ind])))
smCOGsNames <- NULL

for (i in 1:length(smCOGs)){
    tmp <-as.character(10000+smCOGs[i])
    smCOGsNames <- c(smCOGsNames, paste("SMCOG", tmp, sep = ""))
}

### do this and aldo the old version!

names(smCOGs)<- smCOGsNames
rm(smCOGsNames, ind, i)
##### the element (i,j) is the number of colocalizations cosidering the the ith
##### smCOG as fixed
N <- length(smCOGs)
coloc.count <- Matrix(0,N,N)
rownames(coloc.count)<- names(smCOGs)
colnames(coloc.count)<- names(smCOGs)
ind <- which(upper.tri(coloc.count), arr.ind = T)


iterations <- 1:dim(ind)[1]
chuncks <-split(1:dim(ind)[1], ceiling(seq_along(1:dim(ind)[1])/100000))
N <- length(chuncks)
ind.list <- list()
for(jj in 1:N){
    ind.list[[jj]] <- ind[chuncks[[jj]],]
}
save(ind.list, file="indList.Rdata")

rm(iterations,tmp, chuncks, ind,jj, ind.list)

for(ch in 1:N){
    cores <- detectCores()
    cl <- makeCluster(cores-1)
    registerDoParallel(cl)
    load("indList.Rdata")
    
    ind <- ind.list[[ch]]
    rm(ind.list)
    values <-  foreach(i=1:dim(ind)[1], .combine='c') %dopar%{
    res = count.coloc.symmetric(ind[i,1],ind[i,2],smCOGs,clusters)
    res
}
stopCluster(cl)

coloc.count[ind] <- values
rm(values)

save.image("partial_results.Rdata")

cat("\n", (ch/N)*100,"%")
}



save(coloc.count, file = "coloc_count.Rdata")
coloc.count <- as.matrix(coloc.count)
write.csv(coloc.count, "coloc_count.csv")

