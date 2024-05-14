##Exercio 1


pop= read.delim2("1kgp_genomes_metadata.txt.txt", sep= "\t")
AFR= pop[pop$Population == "AFR", ]
AFR_Sample= write.table(AFR, "AFR_SAMPLE.txt", row.names= FALSE, sep= "\t", quote = FALSE)

data= read.table("MEFE_FRQ.frq", header= TRUE, sep = "\t")
counts <- count.fields("MEFE_FRQ.frq", sep = "\t")
table(counts)
novo= merge(data, clinvar, by="Position", all= FALSE)


freq_REF= data$X.ALLELE.REF.FREQ.
freq_ALT= data$X.ALLELE.ALT.FREQ.

#categorizar= function(freq_REF, freq_ALT) {
#  ifelse(freq_REF < 0.01 & freq_ALT < 0.01, "raro",
#        ifelse(freq_REF >= 0.01 | freq_ALT >= 0.01, "comum", NA))
#}


categorizar_freq <- function(freq_ALT) {
  ifelse(freq_ALT < 0.01, "raro",
         ifelse(freq_ALT >= 0.01, "comum", NA))
}


rm(dados)
ALT_freq= categorizar_freq(data$X.ALLELE.ALT.FREQ.)
table(ALT_freq)
data$Classificação= ALT_freq

ALT_frequencia= write.table(data, "frequencia_ALT_MEF2C.txt", sep= "\t", quote= F, row.names= FALSE)

snps= merge(gnomad, data, by= "Position", all = FALSE)
table(snps$Variation.ID)
SNPS_clin= merge(snps, clinvar, by="Variation.ID", all= FALSE)
SNPS_clinP= merge(snps, clinvar, by="Position", all= FALSE)


##PCA
library(tidyverse)
library(ggplot2)

dados_pca <- read.table("mef2cess_pca.eigenvec", header = FALSE)
colnames(dados_pca) <- c("ID", "ID", paste0("PC", 1:10))
dados_pca1=dados_pca[, c(1,3:12)]

# Carregar dados de ancestralidade genômica, se disponíveis
ancestralidade <- read.table("ancestralidade.txt", header = TRUE)


dados_pca_ance <- merge(dados_pca1, ancestralidade, by = "ID", all.x = TRUE)

# Plotar o gráfico de dispersão
ggplot(dados_pca_ance, aes(x = PC1, y = PC2, color= Ancestralidade)) +
  geom_point() +
  labs(x = "Componente Principal 1", y = "Componente Principal 2", title = "PCA POR ANCESTRALIDADE") +
  theme_minimal()

#Adicionar cores para os grupos de ancestralidade (opcional):
#Se você tiver dados de ancestralidade genômica disponíveis, você pode adicionar cores diferentes para cada grupo no gráfico de dispersão. Supondo que você tenha uma coluna "Ancestralidade" em seus dados de ancestralidade, você pode fazer algo assim:

ggplot(dados_pca1, aes(x = PC1, y = PC2, color = Ancestralidade)) +
  geom_point() +
  labs(x = "Componente Principal 1", y = "Componente Principal 2", title = "PCA POR AnCESTRALIDADE") +
  theme_minimal()
