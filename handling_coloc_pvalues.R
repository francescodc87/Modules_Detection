### considering the maximum pvalue in the both directions

load("coloc_pval.Rdata")
ind <- which(upper.tri(coloc.pval), arr.ind=T)
ind2 <- cbind(ind[,2], ind[,1])
values <- cbind(coloc.pval[ind], coloc.pval[ind2])
values <- apply(values, MARGIN = 1, max)
coloc.pval[ind] <- values
coloc.pval[ind2] <- values
### saving it
save(coloc.pval, file = "coloc_pval_MAX.Rdata")
write.csv(coloc.pval, "coloc_pval_MAX.csv")

#### applying BY correction
values <- p.adjust(values, method = "BY")
coloc.pval[ind] <- values
coloc.pval[ind2] <- values

##saving it
save(coloc.pval, file = "coloc_pval_adj.Rdata")
write.csv(coloc.pval, "coloc_pval_adj.csv")

