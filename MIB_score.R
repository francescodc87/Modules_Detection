
"linMap" <- function(x, from, to){
  (x - min(x)) / max(x - min(x)) * (to - from) + from
}



load("ranks_for_MIBscore.Rdata")
load("filtered_modules_2hits.Rdata")


names(weights) <- colnames(Ranks.Used.Metrics)


MIB.score <- apply(Ranks.Used.Metrics,1,weighted.mean, w=weights)

MIB.score  <- linMap(MIB.score, 1,100)

names.modules <- rownames(Metrics.Modules)#
col.name.bup <- c("MIB score", colnames(Metrics.Modules))
Metrics.Modules <- cbind(MIB.score, Metrics.Modules)
rownames(Metrics.Modules) <- names.modules
colnames(Metrics.Modules) <- col.name.bup

save(Modules, Metrics.Modules, file="final_results_modules.Rdata")

