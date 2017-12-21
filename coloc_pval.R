#defining the function
pval.coloc <- function(i,j,smCOGs,clusters,coloc.count){
    ### I have everything
    library(Rmpfr)
    cogA <- smCOGs[i]
    cogB <- smCOGs[j]
    indA <- which(clusters==cogA, arr.ind = T)
    indB <- which(clusters==cogB, arr.ind = T)
    Ntot <- length(which(!is.na(clusters) & clusters!=cogA))
    N1 <- length(which(!is.na(clusters[indA[,1],]) & clusters[indA[,1],]!=cogA))
    N0 <- Ntot-N1
    Btot <- dim(indB)[1]
    B1 <- length(which(!is.na(clusters[indA[,1],]) & clusters[indA[,1],]==cogB))
    B0 <- Btot-B1
    i <- 0:(coloc.count-1)
    out <- sum((chooseZ(N0,(Btot-i))*chooseZ(N1,i)))
    out <- out/chooseZ(Ntot,Btot)
    out <- as.numeric(1-out)
    out
}


#### pvalues coloc


library(parallel)
library(foreach)
library(doParallel)
library(Matrix)
##### p-values for coloc

clusters <- as.matrix(read.csv("clustersCOGSnr_trimmed.csv",
                               stringsAsFactors = F, header = T, row.names = 1 ))


ind <- which(!is.na(clusters) & clusters!= "-")
smCOGs <- sort(unique(as.numeric(clusters[ind])))
smCOGsNames <- NULL

for (i in 1:length(smCOGs)){
    tmp <-as.character(10000+smCOGs[i])
    smCOGsNames <- c(smCOGsNames, paste("SMCOG", tmp, sep = ""))
}

names(smCOGs)<- smCOGsNames
rm(smCOGsNames, ind, i)

load("coloc_count.Rdata")

coloc.count <- as.matrix(coloc.count)
coloc.count[lower.tri(coloc.count)] <- t(coloc.count)[lower.tri(coloc.count)]
diag(coloc.count) <- 0

ind.zeros <- which(coloc.count==0)
ind <- which(coloc.count!=0, arr.ind=T)

N <- length(smCOGs)

coloc.pval <- matrix(0,N,N)
coloc.pval[ind.zeros] <- 1


iterations <- 1:dim(ind)[1]
chuncks <-split(1:dim(ind)[1], ceiling(seq_along(1:dim(ind)[1])/50000))
N <- length(chuncks)
ind.list <- list()
for(jj in 1:N){
    ind.list[[jj]] <- ind[chuncks[[jj]],]
}
save(ind.list, file="indList.Rdata")

rm(iterations,tmp, chuncks, ind,jj, ind.list, ind.zeros)



for(ch in 1:N){
    cores <- detectCores()
    cl <- makeCluster(cores-1)
    registerDoParallel(cl)
    load("indList.Rdata")
    ind <- ind.list[[ch]]
    rm(ind.list)
    values <-  foreach(i=1:dim(ind)[1], .combine='c') %dopar%{
        res = pval.coloc(ind[i,1],ind[i,2],smCOGs,clusters,coloc.count[ind[i,1],ind[i,2]])
        res
    }
    stopCluster(cl)
    
    coloc.pval[ind] <- values
    rm(values)
    save.image("partial_results.Rdata")
    
    cat("\n", (ch/N)*100,"%")
}



save(coloc.pval, file = "coloc_pval.Rdata")
write.csv(coloc.pval, "coloc_pval.csv")


