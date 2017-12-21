### considering the maximum pvalue in the both directions

load("neigh_pval.Rdata")
ind <- which(upper.tri(neigh.pval), arr.ind=T)
ind2 <- cbind(ind[,2], ind[,1])
values <- cbind(neigh.pval[ind], neigh.pval[ind2])
values <- apply(values, MARGIN = 1, max)
neigh.pval[ind] <- values
neigh.pval[ind2] <- values
### saving it
save(neigh.pval, file = "neigh_pval_MAX.Rdata")
write.csv(neigh.pval, "neigh_pval_MAX.csv")

#### applying BY correction
values <- p.adjust(values, method = "BY")
neigh.pval[ind] <- values
neigh.pval[ind2] <- values

##saving it
save(neigh.pval, file = "neigh_pval_adj.Rdata")
write.csv(neigh.pval, "neigh_pval_adj.csv")
