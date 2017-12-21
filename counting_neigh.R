#defining the function
count.neigh <- function(i,j,smCOGs,clusters){
    if(i==j){
        out <- 0
    }else{
        cogA <- smCOGs[i]
        cogB <- smCOGs[j]
        indA <- which(clusters==cogA, arr.ind = T)
        indB <- which(clusters==cogB, arr.ind = T)
        indA.plus <- indA
        indA.minus <- indA
        indA.minus[,2] <-indA[,2]-1
        indA.plus[,2] <-  indA [,2]+1
        indA.minus <- do.call("paste", as.data.frame(indA.minus))
        indA.plus <- do.call("paste", as.data.frame(indA.plus))
        indB <- do.call("paste", as.data.frame(indB))
        out <- length(which(indA.plus%in%indB)) + length(which(indA.minus%in%indB)) 
    }
    out
}



library(parallel)
library(foreach)
library(doParallel)
library(Matrix)
##### Counting neigh interactions

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
##### the element (i,j) is the number of neigh cosidering the the ith
##### smCOG as fixed
N <- length(smCOGs)
neigh.count <- Matrix(0,N,N)
rownames(neigh.count)<- names(smCOGs)
colnames(neigh.count)<- names(smCOGs)
ind <- which(upper.tri(neigh.count), arr.ind = T)


iterations <- 1:dim(ind)[1]
chuncks <-split(1:dim(ind)[1], ceiling(seq_along(1:dim(ind)[1])/100000))
N <- length(chuncks)
ind.list <- list()
for(jj in 1:N){
    ind.list[[jj]] <- ind[chuncks[[jj]],]
}
save(ind.list, file="indList_neigh.Rdata")

rm(iterations,tmp, chuncks, ind,jj, ind.list)

for(ch in 1:N){
    cores <- detectCores()
    cl <- makeCluster(cores-1)
    registerDoParallel(cl)
    load("indList_neigh.Rdata")
    ind <- ind.list[[ch]]
    rm(ind.list)
    values <-  foreach(i=1:dim(ind)[1], .combine='c') %dopar%{
        res = count.neigh(ind[i,1],ind[i,2],smCOGs,clusters)
        res
    }
    stopCluster(cl)
    
    neigh.count[ind] <- values
    rm(values)
    save.image("partial_results_neigh.Rdata")
    
    cat("\n", (ch/N)*100,"%")
}



save(neigh.count, file = "neigh_count.Rdata")
neigh.count <- as.matrix(neigh.count)
write.csv(neigh.count, "neigh_count.csv")

