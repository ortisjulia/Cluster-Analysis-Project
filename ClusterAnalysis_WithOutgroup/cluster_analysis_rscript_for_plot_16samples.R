tbl=read.table("ProjTaxa.1.Q")
par(mar=c(6,4,4,4))
barplot(t(as.matrix(tbl)), names.arg=c("8N05240", "8N05890", "8N06612", "8N73248", "8N73604", "K006", "K010", "K011", "K015", "K019", "Lesina_280", "Lesina_281", "Lesina_282", "Lesina_285", "Lesina_286", "Naxos2"), las=2, main="Structure-like cluster analysis", col="chartreuse3", xlab='Taxa', ylab="Ancestry", border=NA)

