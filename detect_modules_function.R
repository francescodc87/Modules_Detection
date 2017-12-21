library(Matrix)
library(igraph)

detect.modules <- function(c.pval, n.pval, thr, N=44){
  c.pval[c.pval>=thr] <- NA
  c.pval[c.pval<thr] <- 1
  c.pval[is.na(c.pval)] <- 0
  n.pval[n.pval>thr] <- NA
  n.pval[n.pval<=thr] <- 1
  n.pval[is.na(n.pval)] <- 0
  
  bin.matrix <- c.pval+n.pval
  rm(c.pval, n.pval)
  bin.matrix[bin.matrix==2] <- 1
  #bin.matrix <- Matrix(bin.matrix)
  if(any(bin.matrix[upper.tri(bin.matrix)]==0)){
    g1 <- graph_from_adjacency_matrix(adjmatrix = bin.matrix, mode = "undirected",
                                      weighted = NULL, diag = F, add.rownames = "vertex name")
    cliq.g1 <- max_cliques(g1, min=3)
    if (length(cliq.g1)>0){
      Modules <- matrix(NA,length(cliq.g1), N)
      for(i in 1:length(cliq.g1)){
        tmp <- names(cliq.g1[[i]])
        tmp <- sort(tmp)
        Modules[i,1:length(tmp)] <- tmp
        
      }
      Modules<- unique(Modules)      
    }else{
      Modules <- NULL
    }

    
  }else{
    Modules <- NULL
  }
  
    return(Modules)
  
}