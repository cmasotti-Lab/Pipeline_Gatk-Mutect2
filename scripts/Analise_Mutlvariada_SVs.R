
# CARREGANDO DADOS DA TABELA dados GENADOS NO SCRIPT ANTERIOR
# Organizando os dados continuos em (numeric) 
# E os dados categoricos em (factor)
dados0<-dados
#dados0$Response1 <-ifelse(dados0$Response1 == "nCRT-R", 1, 0) # USAR VITAL.STATUS OU RESPOSTA, RESPOSTA2
dados0$Clones<-as.numeric(dados0$Clones)
dados0$Tumor.size<-as.numeric(dados0$Tumor.size)
dados0$Baseline.T<-as.factor(dados0$Baseline.T)
dados0$Baseline.N<-as.factor(dados0$Baseline.N)
dados0$Avoid.TME..1.yes..2.no.<-as.factor(dados0$Avoid.TME..1.yes..2.no.)
#dados1$Avoid.TME..1.yes..2.no. <- gsub("2", "0", dados1$Avoid.TME..1.yes..2.no.)
# Identificar colunas de caractere e transformando em factor
colunas_caractere <- sapply(dados0, is.character)
# Loop para converter colunas de caractere para fator
for (coluna in names(dados0)[colunas_caractere]) {
  dados0[[coluna]] <- as.factor(dados0[[coluna]])
}

# Selecionando principais colunas
dados0<- subset(dados0, 
                select=c("Response1","MATH_score","TMB","MAF_median",
                         "DriverGene","dMMR",
                         "Tratamento","Baseline.T","Baseline.N","Final.Clinical.Stage","Avoid.TME..1.yes..2.no.",
                         "Distance.from.anal.verge","Gender", "Age","Disease.free..1.no..2.yes.","Tumor.size",
                         "Clones","HED_Score","Heterozygosis","Mean_class","HLA.mut",
                         "NAL","nal_bin","NAL_snvs","NAL_indel",
                         "BND","DEL","DUP","INV","SV",
                         "status_BND","status_DEL","status_DUP","status_INV","status_SV"))

# SELECIONANDO VARIAVEIS COM INFORMAÇOES COMPLETAS PARA AS 52 AMOSTRAS
dados1<- subset(dados0, 
                select=c("Response1","MATH_score", "TMB","MAF_median",
                         "DriverGene","dMMR",
                         "Tratamento","Baseline.T","Baseline.N","Gender", "Age",
                         "HED_Score","Heterozygosis","Mean_class","HLA.mut",
                         "NAL", "nal_bin", "NAL_snvs",
                         "BND","DEL","INV",
                         "status_DEL","status_INV","status_SV"))

# Ajuste de um modelo de regressão logística multivariada com todas as variaveis
# Aqui estou testando dois modelos, mas por enquanto eu escolhi o lm()
modelo_logistico <- lm(Response1 ~ .  , data = dados1)
modelo_logistico_generico <- glm(Response1 ~ . , data= dados1, family= binomial)
# Visualizar resumo do modelo
summary(modelo_logistico)
summary(modelo_logistico_generico)
#==============================================================================#
# FOREST PLOT ####
library("forestmodel")
print(forest_model(modelo_logistico))
#==============================================================================#
# Fazendo a regressao liner bivariadas, ou seja, uma variante por vez
summary(lm(Response1 ~ MATH_score, data = dados1))
summary(glm(Response1 ~ MATH_score, data = dados1))

#==============================================================================#
# Encontrando variaveis altamente correlacionadas
#==============================================================================#
library(car)
vif_modelo <- vif(modelo_logistico)
vif_modelo

# refazendo regressão linear multivariada sem as variaveis (TMB, NAL, NAL_snvs)
#==============================================================================#
dados2<- subset(dados1, 
                select=c("Response1","MATH_score","MAF_median",
                         "DriverGene","dMMR",
                         "Tratamento","Baseline.T","Baseline.N","Gender", "Age",
                          "Heterozygosis","Mean_class","HLA.mut",
                         "nal_bin","INV"))

# Ajuste de um modelo de regressão logística multivariada com todas as variaveis
# Aqui estou testando dois modelos, mas por enquanto eu escolhi o lm()
modelo_logistico2 <- lm(Response1 ~ .  , data = dados2)
modelo_logistico_generico2 <- glm(Response1 ~ . , data= dados2, family= binomial)
# Visualizar resumo do modelo
summary(modelo_logistico2)
#summary(modelo_logistico_generico2)

vif(modelo_logistico2)
#vif(modelo_logistico_generico2)
#==============================================================================#
# FOREST PLOT ####
library("forestmodel")
print(forest_model(modelo_logistico2))
#==============================================================================#


#==============================================================================#
# Fazendo a regressao liner bivariadas, ou seja, uma variante por vez
summary(lm(Response1 ~ MATH_score_pROC, data = dados))
summary(lm(Response1 ~ MATH_score, data = dados1))
