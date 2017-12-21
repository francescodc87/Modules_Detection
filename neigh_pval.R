#defining the function
pval.neigh <- function(i,j,smCOGs,clusters,neigh.count){
    ### I have everything
    library(Rmpfr)
    cogA <- smCOGs[i]
    cogB <- smCOGs[j]
    indA <- which(clusters==cogA, arr.ind = T)
    indB <- which(clusters==cogB, arr.ind = T)
    indAplus <- indA
    indAplus[,2] <- indA[,2]+1
    indAminus <- indA
    indAminus[,2] <- indAminus[,2] - 1
    Ntot <- length(which(!is.na(clusters) & clusters!=cogA))
    Atot <- dim(indA)[1]
    Btot <- dim(indB)[1]
    Aedge <- 0
    B1 <- 0
    for(i in 1:Atot){
        if(indA[i,2]==1){
            Aedge <- Aedge+1
            if(clusters[indA[i,1],2] == cogB) {B1 <- B1+1}
        }else{
            if(is.na(clusters[indA[i,1], indA[i,2]+1])){
                Aedge <- Aedge+1
            }else{
                if(clusters[indA[i,1],(indA[i,2]+1)] == cogB) {B1 <- B1+1}
                if(clusters[indA[i,1],(indA[i,2]-1)] == cogB) {B1 <- B1+1}
            }
        }
    }
    N1 <- Aedge+2*(Atot-Aedge)
    Both  <- intersect(indA[,1], indB[,1])
    N0 <- Ntot-N1
    B0 <- Btot-B1
    i <- 0:(neigh.count-1)
    out <- sum((chooseZ(N0,(Btot-i))*chooseZ(N1,i)))
    out <- out/chooseZ(Ntot,Btot)
    out <- as.numeric(1-out)
    out
}



#### pvalues neigh


library(parallel)
library(foreach)
library(doParallel)
library(Matrix)
##### p-values for neigh

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


load("neigh_count.Rdata")

neigh.count <- as.matrix(neigh.count)
neigh.count[lower.tri(neigh.count)] <- t(neigh.count)[lower.tri(neigh.count)]
diag(neigh.count) <- 0

ind.zeros <- which(neigh.count==0)
ind <- which(neigh.count!=0, arr.ind=T)

N <- length(smCOGs)

neigh.pval <- matrix(0,N,N)
neigh.pval[ind.zeros] <- 1

iterations <- 1:dim(ind)[1]
chuncks <-split(1:dim(ind)[1], ceiling(seq_along(1:dim(ind)[1])/5000))
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
        res = pval.neigh(ind[i,1],ind[i,2],smCOGs,clusters,neigh.count[ind[i,1],ind[i,2]])
        res
    }
    stopCluster(cl)
    
    neigh.pval[ind] <- values
    rm(values)
    save.image("partial_results.Rdata")
    
    cat("\n", (ch/N)*100,"%")
}



save(neigh.pval, file = "neigh_pval.Rdata")
write.csv(neigh.pval, "neigh_pval.csv")


