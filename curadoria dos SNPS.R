#CURADORIA DOS DADOS

gnomad= read.delim2("gnomad.txt", sep = "\t")
clinvar2= read.delim2("clinvar_result_all.txt", sep = "\t")

#gnomad2= gnomad2 %>% rename("Variation.ID"= "rsIDs")
clinvar2= clinvar2 %>% rename("Variation.ID"= "dbSNP.ID")
novo_rs2=merge(gnomad, clinvar2, by= "Variation.ID", all=FALSE)

novo_rs2= novo_rs2[, c(1,3,5:6,8:15,19,30)]
novo_rs2$phenotypes = NA
indices <- grep("Intellectual disability", novo_rs2$Condition.s.)
novo_rs2$phenotypes[indices] <- "Intellectual Disability"
novo_rs2$phenotypes[novo_rs2$Condition.s. == "not specified"] = "not specified,not provided"
novo_rs2$phenotypes[novo_rs2$Condition.s. == "not provided"] = "not specified,not provided"
novo_rs2$phenotypes[novo_rs2$Condition.s. == "Inborn genetic diseases"] = "Inborn genetic diseases"

#1filtro: por phenotypo: DI
novo_rs2.DI= novo_rs2[grep("Intellectual disability", novo_rs2$Condition.s.),]

#2filtro por patogenicidade: patogenico, incerto e provavelmente patogenico
novo_rs2.DI_patho= novo_rs2.DI[grep("Pathogenic|Likely pathogenic|Uncertain significance|Pathogenic/Likely pathogenic", novo_rs2.DI$Germline.classification),] #com 13 variaveis

novo_rs2.DI_patho$freq_predominante= apply(novo_rs2.DI_patho[,c(8:15)], 1, function(x)names(x)[which.max(x)])
write.table(novo_rs2.DI_patho, "variantes_atualizadas_interesse.csv", sep = "\t", quote = FALSE, col.names = TRUE)


snps_communs= merge(rs.DI_patho, novo_rs2.DI_patho, by= "Variation.ID", all= FALSE)
snps_communs= snps_communs %>%  rename("Variante"= "Variation.ID")
write.table(snps_communs, "variantes_ comuns_interesse.txt", sep = "\t", quote = FALSE, col.names = TRUE)
