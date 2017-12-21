
"ranks.method" <- function(Metrics.vector){
  Metrics.vector <- as.numeric(Metrics.vector)
  freq.table<- table(Metrics.vector)
  lev <- sort(unique(Metrics.vector))
  tmp <-0
  for(i in 1:length(freq.table)){
    ind <- which(Metrics.vector==lev[i])
    Metrics.vector[ind] <- tmp+ freq.table[i]/2
    tmp <- tmp+freq.table[i]
  }
  return(Metrics.vector)
}


####
load("filtered_modules_2hits.Rdata")

Used.Metrics <- Metrics.Modules[,c(1:3,5:6,8:18)]

Used.Metrics[,"max pvalues"] <- -log(as.numeric(Used.Metrics[,"max pvalues"])) ## this is for the order

ind <- which(Used.Metrics[,"max pvalues"]!="Inf")
ind2 <- which(Used.Metrics[,"max pvalues"]=="Inf")
Used.Metrics[ind2 ,"max pvalues"] <- max(as.numeric(Used.Metrics[ind, "max pvalues"]))+10
rm(ind,ind2)


Ranks.Used.Metrics <-apply(Used.Metrics,2,ranks.method)

save(Ranks.Used.Metrics, file="ranks_for_MIBscore.Rdata")
