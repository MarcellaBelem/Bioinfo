pop=read.delim2("1kgp_genomes_metadata.txt.txt", header= TRUE)
unique(pop$Population)
afr= as.data.frame(pop[pop$Population == "AFR",1])
eas= data.frame(pop[pop$Population == "EAS", 1])
eur= data.frame(pop[pop$Population == "EUR", 1])
sas= data.frame(pop[pop$Population == "SAS", 1])

popu= pop[,1]

write.table(afr, "AFR_sample.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(eas, "EAS_sample.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(eur, "EUR_sample.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(sas, "SAS_sample.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(popu, "all_samples.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

populacao= pop[,c(1,6)]
colnames(populacao)= c("ID", "Ancestralidade")
write.table(populacao, "ancestralidade.txt", sep="\t", quote = FALSE, row.names = FALSE)
