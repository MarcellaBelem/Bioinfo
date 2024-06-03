library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyverse)
#bancos de dados utilizados: SNP-nexus, tabelas do clinvar, gnomad

clinvar= read.delim2("clinvar.txt", sep = "\t")

gnomad= read.table("gnomad.txt", sep = "/t")

sift= read.delim2("sift.txt", sep= "\t")

refseq= read.delim2("refseq.txt", sep = "\t")

rm()
#merge das tabelas clinvar e gnomad pela variantID - o rs

rs= merge(gnomad, clinvar, by= "Variation.ID", all= F)
rs$phenotypes= NA
write.table(rs, "Variantes_MEF2C.txt", sep = "\t", col.names = T, quote = FALSE)
rs= read.table("Variantes_MEF2C.txt", sep = "\t", stringsAsFactors = TRUE)


#1filtro: por phenotypo: DI
rs.DI= rs[grep("Mental retardation|Intellectual Disability", rs$Phenotypes),c(1:15,19:22)]

#2filtro por patogenicidade: patogenico, incerto e provavelmente patogenico
rs.DI_patho= rs.DI[grep("Pathogenic|Likely pathogenic|Uncertain significance", rs.DI$Clinical.Significance),]


write.table(rs.DI_patho, "Variantes_MEF2C_DI.txt", sep= "\t", col.names = T, quote = FALSE)


#grafico representado os achados
rs.DI_patho$Clinical.Significance <- factor(rs.DI_patho$Clinical.Significance, levels = c("Pathogenic", "Likely pathogenic", "Uncertain significance"))
glimpse(rs.DI_patho)
#rs.DI_patho$Clinical.Significance <- factor(rs.DI_patho$Clinical.Significance, levels = c("")

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


rs= read.table("Variantes_MEF2C.txt", sep = "\t", stringsAsFactors = TRUE)
rs$predominant_pop= apply(rs[, c(8:15)], 1, function(x)names(x)[which.max(x)])
tab= rs[rs$phenotypes == "Intellectual Disability",c(1,3,5,6,19:20,22)]
tab= tab[grep("Pathogenic|Likely pathogenic|Uncertain significance", tab$Clinical.Significance),]
unique(rs$Clinical.Significance)
tab2= rs[,c(1,3,5,6,19:20,22,23)]
write.table(tab2, "Variantes_MEF2C_tab.txt", quote= FALSE, sep= "\t", col.names = T )
tab= read.delim2("Variantes_MEF2C_tab.txt", sep = "\t")

glimpse(rs)

frequencias <- table(rs$phenotypes)
frequencias_df <- data.frame(Phenotypes = names(frequencias), Freq = as.vector(frequencias))

# gráfico de pizza: freq fenotipos
library(ggplot2)

ggplot(data = frequencias_df, aes(x = "", y = Freq, fill = Phenotypes)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(round((Freq / sum(Freq)) * 100, 2), "%")), 
            position = position_stack(vjust = 0.5), 
            size = 4, color = "black") +
  theme_void()


frequencias1 <- table(rs$Clinical.Significance)
frequencias1_df <- data.frame(Clinical.Significance = names(frequencias1), Freq = as.vector(frequencias1))

# gráfico de pizza: freq significado variantes
library(ggplot2)

ggplot(data = frequencias1_df, aes(x = "", y = Freq, fill = Clinical.Significance)) +
  geom_bar(stat = "identity", width = 2) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(round((Freq / sum(Freq)) * 100, 2), "%")), 
            size = 4, color = "black", position = position_stack(vjust = 0.5)) +
  theme_void()

##ancestralidade
ancestralidade= rs[, c(1,8:15,20,22)]
ancestralidade= ancestralidade %>% rename( Variante= Variation.ID, AFR = AFR.Frequency,
       AMR = AMR.Frequency,
       ASJ = ASJ.Frequency,
       EAS = EAS.Frequency,
       FIN = FIN.Frequency,
       NFE = NFE.Frequency,
       OTH = OTH.Frequency,
       SAS = SAS.Frequency, 
       Clinical_snp= Clinical.Significance, 
       Phenotypes = phenotypes) %>% 
  mutate(Doença = if_else(Phenotypes == "Intellectual Disability", as.factor("YES"), as.factor("NO")))
ancestralidade$freq_predominante= apply(ancestralidade[,-1], 1, function(x)names(x)[which.max(x)])

#só comas as 5 variantes
ancestralidadeDI=  novo_rs2[, c(1,5:12, 14,15)] 
ancestralidadeDI= merge(ancestralidadeDI, snps_communs, by= "Variation.ID", all=FALSE)
ancestralidadeDI= ancestralidadeDI[,c(1:9,11,27)]%>% rename( Variante= Variation.ID, AFR = AFR.Frequency,
                                           AMR = AMR.Frequency,
                                           ASJ = ASJ.Frequency,
                                           EAS = EAS.Frequency,
                                           FIN = FIN.Frequency,
                                           NFE = NFE.Frequency,
                                           OTH = OTH.Frequency,
                                           SAS = SAS.Frequency, 
                                           Clinical_snp= Clinical.Significance,
                                           Phenotypes = phenotypes) 
ancestralidadeDI$freq_predominante= apply(ancestralidadeDI[,-1], 1, function(x)names(x)[which.max(x)])
ancestralidadeDI= ancestralidadeDI[ancestralidadeDI$Phenotypes == "Intellectual Disability", ] 
ancestralidadeDI= ancestralidadeDI[ancestralidadeDI$Clinical_snp == "Uncertain significance" | ancestralidadeDI$Clinical_snp == "Pathogenic", ]
ancestralidadeDI$Clinical_snp <- factor(ancestralidadeDI$Clinical_snp, levels = c("Pathogenic", "Likely pathogenic", "Uncertain significance","Likely benign", "Benign", "Conflicting interpretations of pathogenicity"))


ggplot(ancestralidadeDI, aes(x = Variante, y = freq_predominante, fill = freq_predominante)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_jitter(aes(x = Variante, y = " ", color = Clinical_snp, shape = Clinical_snp), position = position_jitter(width = 0.2), size= 2) +
  scale_color_manual(values = c("red","blue", "green", "purple", "orange")) +
  labs(fill = "Population",
       x = "Variants", y = "Frequencies", shape= "Clinical Significance", color= "Clinical Significance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

#separar de 20 a 20 variantes
unique(ancestralidade$Phenotypes)
cores= c("Pathogenic"="red3", "Likely benign"= "steelblue","Benign"=  "green2","Uncertain significance"="purple1", "Conflicting interpretations of pathogenicity"= "orange2","Likely pathogenic"= "pink1",
         "Intellectual Disability"= "red", "not specified,not provided"= "blue", "History of neurodevelopmental disorder"= "green")
ancestralidade1= ancestralidade[1:20,]
ggplot(ancestralidade1, aes(x = Variante, y = freq_predominante, fill = Clinical_snp)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_jitter(aes(x = Variante, y = " ", color = Phenotypes, shape = Phenotypes), position = position_jitter(width = 0.3), size= 2) +
  scale_fill_manual(values = cores) +
  scale_color_manual(values = cores) +
  labs(fill = "Clinical Significance",
       x = "Variants", y = "Frequencies") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ancestralidade2= ancestralidade[21:40,]
ggplot(ancestralidade2, aes(x = Variante, y = freq_predominante, fill = Clinical_snp)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_jitter(aes(x = Variante, y = " ", color = Phenotypes, shape = Phenotypes), position = position_jitter(width = 0.4), size= 2) +
  scale_fill_manual(values = cores) +
  scale_color_manual(values = cores) +
  labs(fill = "Clinical Significance",
       x = "Variants", y = "Frequencies") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ancestralidade3= ancestralidade[41:61,]
ggplot(ancestralidade3, aes(x = Variante, y = freq_predominante, fill = Clinical_snp)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_jitter(aes(x = Variante, y = " ", color = Phenotypes, shape = Phenotypes), position = position_jitter(width = 0.3), size= 2) +
  scale_fill_manual(values = cores) +
  scale_color_manual(values = cores) +
  labs(fill = "Clinical Significance",
       x = "Variants", y = "Frequencies") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#################
#associação com o gwas
gwas= read.csv("gwas_mef2c.txt", sep = "\t")
vcf=read.delim2("MEFE_AllFRQ.frq", sep = "\t")
common_snps <- intersect(gwas$Position, vcf$Position)
association= merge(gwas, vcf, by= "Position", all= FALSE)
common_snps <- intersect(association$Variation.ID, clinvar$Variation.ID)

# testes estatisticos
#
frequencias <- table(rs$phenotypes)
frequencias_df <- data.frame(Phenotypes = names(frequencias), Freq = as.vector(frequencias))
chisq.test(frequencias_df$Freq)

#
frequencias1 <- table(rs$Clinical.Significance)
frequencias1_df <- data.frame(Clinical.Significance = names(frequencias1), Freq = as.vector(frequencias1))
chisq.test(frequencias1_df$Freq)


##teste todas as amostras
#com todas as variantes
dados_test1= ancestralidade[,1:9]
dados_longos1 <- reshape2::melt(dados_test1, id.vars = "Variante", variable.name = "Populacao", value.name = "Frequencia")
dados_longos1$Frequencia <- as.numeric(dados_longos1$Frequencia)
dados_longos1$Frequencia <- dados_longos1$Frequencia * 100
shapiro.test(dados_longos1$Frequencia) #anormal

#kruskalis-test
kruskal_test_result1 <- kruskal.test(Frequencia ~ Populacao, data = dados_longos1)
#Há uma diferença significativa entre as populações

ggplot(dados_longos1, aes(x = Populacao, y = Frequencia, fill = Populacao)) +
  geom_bar(stat = "identity") +
  labs(x = "Population",
       y = "Frequency") +
  theme_minimal()

freq_clinico=xtabs(~Clinical_snp+freq_predominante, ancestralidade)
resultado_teste <- chisq.test(freq_clinico) #não significativo
resultados_testeFisher = fisher.test(freq_clinico)

barplot(freq_clinico, col = 2:6, beside = TRUE, legend.text = TRUE,
        main = "",
        xlab = "Frequency Predominante", ylab = "Frequency",
        args.legend = list(x = "topright", y = max(freq_clinico) + 100, bty = "n", cex = 0.5))
       
###teste t para ver se havia diferença significativa entre as populações entre as 5 variantes
dados_test= ancestralidadeDI[,1:9]
str(dados_test)
dados_test[,-1] <- apply(dados_test[,-1], 2, as.numeric)
dados_test[,-1] <- lapply(dados_test[,-1], function(x) x * 100)


# Calcular as diferenças intergênicas usando um teste t de Student simples como exemplo
  resultados_teste <- apply(dados_test[, -1], 1, function(x) t.test(x)$p.value)
  
# Definir um limiar de significância
  limiar <- 0.05
  
# Identificar variantes com diferenças significativas
variantes_significativas <- which(resultados_teste < limiar)
#não há diferenças significativas entre as frequências de cada variantes nas diferentes populações

resultados_df <- data.frame(
  Variante = dados_test$Variante,
  P_value = resultados_teste,
  Significativa = resultados_teste < limiar
)

# Plotar os resultados usando ggplot2
library(ggplot2)
ggplot(resultados_df, aes(x = Variante, y = -log10(P_value), fill = Significativa)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
  labs(title = "Significância das Variantes",
       x = "Variante",
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

 
library(tidyr)
library()

dados_longos <- reshape2::melt(dados_test, id.vars = "Variante", variable.name = "Populacao", value.name = "Frequencia")
dados_longos$Frequencia <- as.numeric(dados_longos$Frequencia)
dados_longos$Frequencia <- dados_longos$Frequencia * 100
shapiro.test(dados_longos$Frequencia) #anormal

#kruskalis-test
kruskal_test_result <- kruskal.test(Frequencia ~ Populacao, data = dados_longos)
#Não há uma diferença significativa entre as populações

ggplot(dados_longos, aes(x = Populacao, y = Frequencia, fill = Populacao)) +
  geom_bar(stat = "identity") +
  labs(x = "Population",
       y = "Frequency") +
  theme_minimal()


#fst
install.packages("hierfstat")
install.packages("adegenet")
library(hierfstat)
library(adegenet)

data= ancestralidadeDI[,1:9]
freq_data= data[ ,-1]

calculate_fst <- function(freq_data) {
  fst_values <- numeric(nrow(freq_data))
  
  for (i in 1:nrow(freq_data)) {
    freqs <- as.numeric(freq_data[i, ])
    
    # Média das frequências alélicas entre populações
    p_bar <- mean(freqs)
    
    # Heterozigosidade esperada dentro das populações (HS)
    hs <- mean(2 * freqs * (1 - freqs))
    
    # Heterozigosidade esperada total (HT)
    ht <- 2 * p_bar * (1 - p_bar)
    
    # FST
    if (ht != 0) {
      fst_values[i] <- (ht - hs) / ht
    } else {
      fst_values[i] <- 0
    }
  }
  
  return(fst_values)
}

# Calcular FST
data$FST <- calculate_fst(freq_data)
print(data)
#apresenta alta diferenciação genetica entre as populações

data_long <- reshape2::melt(data, id.vars = "Variante")
colnames(data_long) <- c("Variante", "Populacao", "Frequencia")
data_long$Frequencia= as.numeric(data_long$Frequencia)

ggplot(data_long, aes(x = Variante, y = Frequencia, fill = Populacao)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(Frequencia, 2)), vjust = -0.3, position = position_dodge(0.9), size = 3.5) +
  labs(title = "Frequências Alélicas por Variante e População",
       x = "Variante",
       y = "Frequência Alélica") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data, aes(x = Variante, y = FST, fill = Variante)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(FST, 2)), vjust = -0.3, size = 3.5) +
  labs(title = "Valores de FST por Variante",
       x = "Variante",
       y = "FST") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

