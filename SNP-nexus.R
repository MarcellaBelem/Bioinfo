library(tidyverse)
#bancos de dados utilizados: SNP-nexus, tabelas do clinvar, gnomad

clinvar= read.delim2("clinvar.txt", sep = "\t")

gnomad= read.delim2("gnomad.txt", sep = "\t")

sift= read.delim2("sift.txt", sep= "\t")

refseq= read.delim2("refseq.txt", sep = "\t")

rm()
#merge das tabelas clinvar e gnomad pela variantID - o rs

rs= merge(gnomad, clinvar, by= "Variation.ID", all= F)
rs$phenotypes= NA
write.table(rs, "Variantes_MEF2C.txt", sep = "\t", col.names = T, quote = FALSE)
rss= read.table("Variantes_MEF2C.txt", sep = "\t")

#1filtro: por phenotypo: DI
rs.DI= rs[grep("Mental retardation|Intellectual Disability", rs$Phenotypes),c(1:15,19:22)]

#2filtro por patogenicidade: patogenico, incerto e provavelmente patogenico
rs.DI_patho= rs.DI[grep("Pathogenic|Likely pathogenic|Uncertain significance", rs.DI$Clinical.Significance),]


#grafico representado os achados
rs.DI_patho$Clinical.Significance <- factor(rs.DI_patho$Clinical.Significance, levels = c("Pathogenic", "Likely pathogenic", "Uncertain significance"))
glimpse(rs.DI_patho)
#rs.DI_patho$Clinical.Significance <- factor(rs.DI_patho$Clinical.Significance, levels = c("")

#1 por variante tipo

pie(table(rs.DI_patho$Clinical.Significance), col= rainbow(length(levels(rs.DI_patho$Clinical.Significance)))) 
legend("bottomright", levels(rs.DI_patho$Clinical.Significance), title = "Significado Clínico", fill = rainbow(length(levels(rs.DI_patho$Clinical.Significance))), cex= 0.8)

tabela0= table(rs.DI_patho$Clinical.Significance)
porcentagem0=  round(prop.table(tabela0) * 100, 2)
pie(porcentagem0, labels = paste(names(porcentagem0), "(", porcentagem0, "%)", sep = ""), col= rainbow(length(levels(rs.DI_patho$Clinical.Significance))), cex= 0.8)
pie(porcentagem0, labels = paste(names(porcentagem0), "(", porcentagem0, "%)", sep = ""), col= rainbow(length(levels(rs.DI_patho$phenotypes))), cex= 0.8)
#ver no ggplot
tabelavari= data.frame(xtabs(~Clinical.Significance+phenotypes, rs))

 ggplot(tabela0, x= "", y= phenotypes, fil= Clinical.Significance)
 geom_bar(stat = "identity", width = 1) +
   coord_polar("y", start = 0) +
   geom_text(aes(label = paste0(round(percentagem), "%")), position = position_stack(vjust = 0.5)) +
   theme_void() +
   labs(title = "Gráfico de Pizza com Porcentagens")

#2 por phenotype
rs.DI_patho <- factor(rs.DI_patho$phenotypes, levels = c("Intellectual Disability", "Mental retardation")
tabela2= table(rs.DI_patho$phenotypes)
  porcentagem2=  round(prop.table(tabela2) * 100, 2)
pie(porcentagem0, labels = paste(names(porcentagem0), "(", porcentagem0, "%)", sep = ""), col= rainbow(length(levels(rs.DI_patho$phenotypes))), cex= 0.8)
legend("top", levels(rs.DI_patho$Clinical.Significance), title = "DI vs significado clinico da variante", fill = rainbow(length(levels(rs.DI_patho$phenotypes))), cex= 0.8)

#barplot( table(rs.DI_patho$Clinical.Significance), col= rainbow(length(levels(rs.DI_patho$Clinical.Significance))))
#legend("top", levels(rs.DI_patho$Clinical.Significance), title = "Significado Clínico", fill = rainbow(length(levels(rs.DI_patho$Clinical.Significance))), cex= 0.8)
                      
#curadoria dos resultados no clinvar; escolher 3 variantes e ver se os resultados correspondem em outro site
                      
                      
unique(rs.DI$Clinical.Significance)
#grafico total snp por phenotypo e clinica variante
rs$phenotypes[rs$Phenotypes == "not specified,Intellectual Disability, Stereotypic Movements, Epilepsy, and/or Cerebral Malformations,not provided"] = "Intellectual Disability"
rs$phenotypes[rs$Phenotypes == "not specified,Intellectual Disability, Stereotypic Movements, Epilepsy, and/or Cerebral Malformations"] = "Intellectual Disability" 
rs$phenotypes[rs$Phenotypes == "not specified,not provided"] = "not specified,not provided"
rs$phenotypes[rs$Phenotypes == "not specified"] = "not specified,not provided"
rs$phenotypes[rs$Phenotypes == "not provided"] = "not specified,not provided"    
rs$phenotypes[rs$Phenotypes == "History of neurodevelopmental disorder,not specified,not provided"] = "History of neurodevelopmental disorder"
rs$phenotypes[rs$Phenotypes == "History of neurodevelopmental disorder,not specified"] = "History of neurodevelopmental disorder"
rs$phenotypes[rs$Phenotypes == "Intellectual Disability, Stereotypic Movements, Epilepsy, and/or Cerebral Malformations"] = "Intellectual Disability"
rs$phenotypes[rs$Phenotypes == "Mental retardation, stereotypic movements, epilepsy, and/or cerebral malformations,not provided"] = "Intellectual Disability"
rs$phenotypes[rs$Phenotypes == "Mental retardation, stereotypic movements, epilepsy, and/or cerebral malformations"] = "Intellectual Disability"
rs$phenotypes[rs$Phenotypes == "Intellectual Disability, Stereotypic Movements, Epilepsy, and/or Cerebral Malformations"] = "Intellectual Disability"
rs$phenotypes[rs$Phenotypes == "Mental retardation, stereotypic movements, epilepsy, and/or cerebral malformations,not specified,not provided"] = "Intellectual Disability"


rs$phenotypes= factor(rs$phenotypes, levels = c("Intellectual Disability", "History of neurodevelopmental disorder", "not specified,not provided"))

rs$Clinical.Significance[rs$Clinical.Significance == "Benign/Likely benign"] = "Likely benign"
rs$Clinical.Significance <- factor(rs$Clinical.Significance, levels = c("Pathogenic", "Likely pathogenic", "Uncertain significance","Likely benign", "Benign", "Conflicting interpretations of pathogenicity"))
unique(rs$Clinical.Significance)

##legend("bottom", levels(rs$phenotypes), title = "Fenótipo", fill = rainbow(length(levels(rs$phenotypes))), cex= 0.6)

##variancia clinica
tabela= table(rs$Clinical.Significance)
porcentagem=  round(prop.table(tabela) * 100, 2)
pie(porcentagem, labels = paste(names(porcentagem), "(", porcentagem, "%)", sep = ""), col= rainbow(length(levels(rs$Clinical.Significance))), cex= 0.7)
legend("topright", legend = levels(rs$Clinical.Significance), fill = rainbow(length(levels(rs$Clinical.Significance))), cex = 0.8, title = "Clinical Significance")

vari= data.frame(tabela, porcentagem)

#phenotypos
tabela1= table(rs$phenotypes)
porcentagem1=  round(prop.table(tabela1) * 100, 2)
pie(porcentagem1, labels = paste(names(porcentagem1), "(", porcentagem1, "%)", sep = ""), col= rainbow(length(levels(rs$phenotypes))), cex= 0.8)

pheno= data.frame(tabela1, porcentagem1)
unique(rs$Phenotypes)




## frequencias
data= rs[,c(1,8:15)]
data$AFR.Frequency= as.numeric(data$AFR.Frequency)
data$AMR.Frequency= as.numeric(data$AMR.Frequency )
data$ASJ.Frequency= as.numeric(data$ASJ.Frequency)
data$EAS.Frequency= as.numeric(data$EAS.Frequency)
data$FIN.Frequency= as.numeric(data$FIN.Frequency)
data$NFE.Frequency= as.numeric(data$NFE.Frequency)
data$OTH.Frequency= as.numeric(data$OTH.Frequency)
data$SAS.Frequency= as.numeric(data$SAS.Frequency)


intervalos= cut(data$SAS.Frequency, breaks = c(-Inf, 0.01, Inf),labels= c("<0.01", ">0.01"))
data$SAS.Frequency= factor(intervalos, levels= c("<0.01", ">0.01"))

intervalos= cut(data$OTH.Frequency, breaks = c(-Inf, 0.01, Inf),labels= c("<0.01", ">0.01"))
data$OTH.Frequency= factor(intervalos, levels= c("<0.01", ">0.01"))

intervalos= cut(data$EAS.Frequency, breaks = c(-Inf, 0.01, Inf),labels= c("<0.01", ">0.01"))
data$EAS.Frequency= factor(intervalos, levels= c("<0.01", ">0.01"))

intervalos= cut(data$ASJ.Frequency, breaks = c(-Inf, 0.01, Inf),labels= c("<0.01", ">0.01"))
data$ASJ.Frequency= factor(intervalos, levels= c("<0.01", ">0.01"))

intervalos= cut(data$NFE.Frequency, breaks = c(-Inf, 0.01, Inf),labels= c("<0.01", ">0.01"))
data$NFE.Frequency= factor(intervalos, levels= c("<0.01", ">0.01"))

intervalos= cut(data$FIN.Frequency, breaks = c(-Inf, 0.01, Inf),labels= c("<0.01", ">0.01"))
data$FIN.Frequency= factor(intervalos, levels= c("<0.01", ">0.01"))

intervalos= cut(data$AFR.Frequency, breaks = c(-Inf, 0.01, Inf),labels= c("<0.01", ">0.01"))
data$AFR.Frequency= factor(intervalos, levels= c("<0.01", ">0.01"))

intervalos= cut(data$AMR.Frequency, breaks = c(-Inf, 0.01, Inf),labels= c("<0.01", ">0.01"))
data$AMR.Frequency= factor(intervalos, levels= c("<0.01", ">0.01"))


freq= xtabs(~Variation.ID+AFR.Frequency, data=data)
freq2= xtabs(~Variation.ID+AMR.Frequency, data=data)
freq3= xtabs(~Variation.ID+EAS.Frequency, data=data)
freq4= xtabs(~Variation.ID+ASJ.Frequency, data=data)
fre5 = xtabs(~Variation.ID+OTH.Frequency, data=data)
freq6= xtabs(~Variation.ID+SAS.Frequency, data=data)
freq7=xtabs(~Variation.ID+FIN.Frequency, data=data)
freq8= xtabs(~Variation.ID+NFE.Frequency, data=data)

freqt= merge(freq,freq2, by= "Variation.ID", all= FALSE)
freqto= merge(freqt,freq3, by= "Variation.ID", all= FALSE)
freqtot= merge(freqto, freq4, by ="Variation.ID", all=FALSE)
freqq= merge(freqtot, fre5, by= "Variation.ID", all= FALSE)
freqq1= merge(freqq, freq6, by= "Variation.ID", all= FALSE)
freqq2= merge(freqq1, freq7, by= "Variation.ID", all= FALSE)
freq_total= merge(freqq2, freq8, by= "Variation.ID", all= FALSE)
freq_total1= freq_total[,c(1, 2,4,6,8,10,12,14,16)]

hist()