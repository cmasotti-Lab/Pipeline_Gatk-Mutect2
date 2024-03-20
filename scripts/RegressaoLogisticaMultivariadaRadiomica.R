library("forestmodel")
library(car)
library(forestplot)
library(data.table)
library("survminer")

# CARREGANDO DADOS DA TABELA dados GENADOS NO SCRIPT ANTERIOR
# Organizando os dados continuos em (numeric) 
# E os dados categoricos em (factor)
# dados0<-dados
# 
# dados0$CEA.inicial<-as.numeric(dados0$CEA.inicial)

# dir.create("D:/PROJETOS-HSL_BP/RegressaoLogisticaMultivariada/")
# fwrite(dados0,"D:/PROJETOS-HSL_BP/RegressaoLogisticaMultivariada/dados0.tsv", quote = F, sep="\t")

dados0 <- fread("D:/PROJETOS-HSL_BP/RegressaoLogisticaMultivariada/dados0.tsv", quote = F, sep="\t")
dados0 <- as.data.frame(lapply(dados0, as.factor))
#==============================================================================#
# Estamos avaliando os pacientes que RESPONDEM ao tratamento ####
#==============================================================================#
dados0$ClinicalOutcome <-as.factor(ifelse(dados0$Response1 == "nCRT-NR", 0, 1)) 

#==============================================================================#
# Selecionando principais colunas
dados0 <- subset(dados0, 
                 select=c("ID_Exoma","ClinicalOutcome","MATH_score", "Treatment.Protocol","Baseline.T",
                          "CEA.inicial", "Baseline.N","Final.Clinical.Stage","Avoid.TME..1.yes..2.no.",
                          "Distance.from.anal.verge","Gender", "Age", "Tumor.size","MATH_score."))

dados0$Tumor.size <- as.numeric(as.character(dados0$Tumor.size))
dados0$MATH_score <- as.numeric(as.character(dados0$MATH_score))
dados0$CEA.inicial <- as.numeric(as.character(dados0$CEA.inicial))
dados0$Distance.from.anal.verge <- as.numeric(as.character(dados0$Distance.from.anal.verge))
dados0$Age<-as.numeric(dados0$Age)

dados0$Baseline.T. <-as.factor(ifelse(dados0$Baseline.T == "2", "2", "3-4"))
dados0$Final.Clinical.Stage. <-as.factor(ifelse(dados0$Final.Clinical.Stage == "3", "3", "1-2"))


#==============================================================================#
# DISPERSSÃO DAS MATH_SCORE DAS AMOSTRAS (R/NR) ####
#==============================================================================#
ggplot(dados0, aes(x=MATH_score, y=ClinicalOutcome)) + 
  geom_point() + 
  stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)

#==============================================================================#
# REGRESSÃO LOGISTICA SIMPLES   ####
#==============================================================================#
modelo_logistico<-glm(ClinicalOutcome ~ MATH_score., family= binomial(link = 'logit'), data = dados0)
modelo_logistico
forest_model(modelo_logistico)+
  theme(text = element_text(size = 14))

#==============================================================================#
variaveis <- c("MATH_score", "Treatment.Protocol", "Baseline.T","Baseline.T.",
               "Baseline.N", "Gender", "Age", "Distance.from.anal.verge", "Tumor.size",
               "Final.Clinical.Stage","Final.Clinical.Stage.",
               "CEA.inicial", "MATH_score.")

results <- data.frame(Variable = character(), P_value = numeric(), stringsAsFactors = FALSE)

for (variable in variaveis) {
  # variable<-"MATH_score."
  cat("=======================================================")
  formula <- as.formula(paste("ClinicalOutcome ~", variable))
  modelo_simples <- glm(formula , family= binomial(link = 'logit'), data = dados0)
  # Extract the p-value from the summary
  p_value <- summary(modelo_simples)$coefficients[, "Pr(>|z|)"][2]
  print(summary(modelo_simples))
  # Store results in a data frame
  results <- rbind(results, data.frame(Variable = variable, P_value = p_value))
  print(forest_model(modelo_simples)+
          theme(text = element_text(size = 14)))
}
# Print the formatted table
knitr::kable(results)

#==============================================================================#
# REGRESSÃO LOGISTICA MULTIVARIADAS ####
#==============================================================================#
# MODELO 1: 52 AMOSTRAS
modelo1<- subset(dados0, 
                 select=c("ClinicalOutcome","MATH_score", 
                          "Treatment.Protocol","Baseline.T.","Baseline.N","Gender", "Age"))
# modelo1$Baseline.T. <-as.factor(ifelse(modelo1$Baseline.T == "2", "2", "3-4"))

# Ajuste de um modelo de regressão logística multivariada com todas as variaveis
modelo_logistico <- glm(ClinicalOutcome ~ .  ,family  = binomial(link = 'logit'), data = modelo1)
summary(modelo_logistico)

# FOREST PLOT ####
print(forest_model(modelo_logistico)+
        theme(text = element_text(size = 14)))

# Encontrando variaveis altamente correlacionadas
vif(modelo_logistico)

#==============================================================================#
# CENARIO 2
modelo2<- subset(dados0, 
                 select=c("ClinicalOutcome","MATH_score", 
                          "Treatment.Protocol","Baseline.T.","Baseline.N","Gender", "Age",
                          "Final.Clinical.Stage","Distance.from.anal.verge","Tumor.size"))
modelo_logistico <- glm(ClinicalOutcome ~ . , family = binomial(link = 'logit'), data = modelo2)
summary(modelo_logistico)
table(modelo2$Final.Clinical.Stage)
modelo2$Final.Clinical.Stage <-as.factor(ifelse(modelo2$Final.Clinical.Stage == "3", "3", "1-2"))
modelo_logistico <- glm(ClinicalOutcome ~ . , family = binomial(link = 'logit'), data = modelo2)
summary(modelo_logistico)
# FOREST PLOT ####
print(forest_model(modelo_logistico)+   
        theme(text = element_text(size = 14)))

modelo2<- subset(modelo2, 
                 select=c("ClinicalOutcome","MATH_score", 
                          "Treatment.Protocol","Baseline.T.","Baseline.N","Gender", "Age",
                          "Distance.from.anal.verge","Tumor.size"))
# Ajuste de um modelo de regressão logística multivariada com todas as variaveis
modelo_logistico <- glm(ClinicalOutcome ~ . , family = binomial(link = 'logit'), data = modelo2)
summary(modelo_logistico)
# FOREST PLOT ####
print(forest_model(modelo_logistico)+ 
        theme(text = element_text(size = 14)))
# Encontrando variaveis altamente correlacionadas
vif(modelo_logistico)


#==============================================================================#
# SELECIONANDO VARIAVEIS COM INFORMAÇOES COMPLETAS PARA AS 52 AMOSTRAS
modelo3<- subset(dados0, 
                 select=c("ClinicalOutcome","MATH_score.", 
                          "Treatment.Protocol","Baseline.T","Baseline.N","Gender", "Age",
                          "Distance.from.anal.verge","Tumor.size"))
modelo3$Baseline.T <-as.factor(ifelse(modelo3$Baseline.T == "2", "2", "3-4"))
# modelo3$Final.Clinical.Stage <-as.factor(ifelse(modelo3$Final.Clinical.Stage == "3", "3", "1-2"))
modelo_logistico <- glm(ClinicalOutcome ~ . , family = binomial(link = 'logit'), data = modelo3)
summary(modelo_logistico)
# FOREST PLOT ####
print(forest_model(modelo_logistico)+   
        theme(text = element_text(size = 14)))

# Encontrando variaveis altamente correlacionadas
vif(modelo_logistico)

#==============================================================================#
# Fazendo a regressao logistica simples do MATH_score classificado
modelo_logistico<-glm(ClinicalOutcome ~ MATH_score. , family = binomial(link = 'logit'), data = dados0)
summary(modelo_logistico)
# FOREST PLOT ####
print(forest_model(modelo_logistico)+
        theme(text = element_text(size = 14)))

#==============================================================================#
col.vars<- "MATH_score."
for (col.var in col.vars) {
  resp <- "ClinicalOutcome"
  # col.var<- "dMMR"
  dt3<- matrix(table(subset(dados0, select = c(col.var, resp))), nr=2)
  if(ncol(dt3)>1){
    resultado<-fisher.test(dt3)
    valor_p <- resultado$p.value
    if(valor_p <= 0.05){
      print("=============================================================")
      print(col.var)
      print(resp)
      print(as.matrix(table(subset(dados0, select = c(col.var, resp))), nr=2))
      print(resultado)
      print("=============================================================")
      print(chisq.test(dt3))
      print("=============================================================")
    }
  }
}
#==============================================================================#






#==============================================================================#
                            # DADOS DE RADIOMICA ####
#==============================================================================#



#==============================================================================#
# CARREGANDO TABELA ENTROPIA, ASSIMETRIA E RADIOMICA ####
#==============================================================================#
DWI.ADC <- read.xlsx("D:/PROJETOS-HSL_BP/Planilha radiomica edit.xlsx", sheetIndex = "DWI e ADC", header = T)
radiomic <- read.xlsx("D:/PROJETOS-HSL_BP/Planilha radiomica edit.xlsx", sheetIndex = "Radiômica", header = T)

Clinical <- read.xlsx("D:/PROJETOS-HSL_BP/Dados Clinicos 2024 ALLversions WXS-TCRseq.xlsx", sheetIndex = "73samples", header = T)
Clinical$ClinicalOutcome <-as.factor(ifelse(Clinical$Response1 == "nCRT-R", 1, 0) )

#==============================================================================#
# EXCLUIR COLUNA PACIENTES e SEM VALORES DE DIFUSAO ####
#==============================================================================#
DWI.ADC<- DWI.ADC[,-c(1,3,14)]
radiomic<- radiomic[,-c(1,3)]
DWI.ADC<- DWI.ADC[!is.na(DWI.ADC$Difusao.Minimum),]
DWI.ADC$Difusão.Média <- as.numeric(as.character(DWI.ADC$Difusão.Média))
radiomic$Volume..mm3. <- as.numeric(as.character(radiomic$Volume..mm3.))

#==============================================================================#
# JUNTANDO TABELAS DADOS CLINICOS COM DWI.ADC ####
#==============================================================================#
DWI.ADC.0 <-merge(subset(Clinical,Response3 != "Metastatic", select = c(2,32)), DWI.ADC, by= "ID_Exoma", all.y = T)
DWI.ADC.0[DWI.ADC.0$ID_Exoma=="ROP-111", ]$ClinicalOutcome <- 0

#==============================================================================#
# JUNTANDO TABELAS DADOS CLINICOS COM RADIOMICA ####
#==============================================================================#
#radiomic.0 <-merge(subset(dados0, select = c(1,2)), radiomic, by= "ID_Exoma")
radiomic.0 <-merge(subset(Clinical,Response3 != "Metastatic", select = c(2,32)), radiomic, by= "ID_Exoma", all.y = T)
radiomic.0[radiomic.0$ID_Exoma=="ROP-111", ]$ClinicalOutcome <- 0


#==============================================================================#
# SELECIONANDO VARIAVEIS
variaveis.DWI.ADC.0 <- colnames(DWI.ADC.0[,-c(1,2)])
variaveis.radiomic.0 <- colnames(radiomic.0[,-c(1,2)])
#==============================================================================#
# SELECIONE O DATASET
variaveis <- variaveis.DWI.ADC.0
dt <- DWI.ADC.0
variaveis <- variaveis.radiomic.0
dt <- radiomic.0
#==============================================================================#

results <- data.frame(Variable = character(), P_value = numeric(), stringsAsFactors = FALSE)

for (variable in variaveis) {
  # variable<-"MATH_score."
  cat("=======================================================")
  cat(variable)
  formula <- as.formula(paste("ClinicalOutcome ~", variable))
  modelo_simples <- glm(formula , family= binomial(link = 'logit'), data = dt)
  # Extract the p-value from the summary
  p_value <- summary(modelo_simples)$coefficients[, "Pr(>|z|)"][2]
  print(summary(modelo_simples))
  # Store results in a data frame
  results <- rbind(results, data.frame(Variable = variable, P_value = p_value))
  print(forest_model(modelo_simples)+
          theme(text = element_text(size = 14)))
}
# Print the formatted table
knitr::kable(results)

#==============================================================================#
# REGRESSÃO LOGISTICA MULTIVARIADAS  ---     DWI.ADC.0   ---####
#==============================================================================#
# MODELO DWI.ADC.0 :22 amostras 
modelo1<- DWI.ADC.0[,-c(1)]
# MODELO radiomic.0 :25 amostras 
modelo1<- radiomic.0[,c(2,5,8,9)]

# Ajuste de um modelo de regressão logística multivariada com todas as variaveis
modelo_logistico <- glm(ClinicalOutcome ~ .  ,family  = binomial(link = 'logit'), data = modelo1)
summary(modelo_logistico)

# FOREST PLOT ####
print(forest_model(modelo_logistico)+
        theme(text = element_text(size = 14)))

# Encontrando variaveis altamente correlacionadas
vif(modelo_logistico)



#==============================================================================#
# REGRESSÃO LOGISTICA MULTIVARIADAS DWI.ADC.0  + Radiomic + Clincal  ---####
#==============================================================================#
# MODELO DWI.ADC.0 :22 amostras 
modelo2<- merge(DWI.ADC.0[,c(1,2,5,8,9)], radiomic.0[,c(1,5,8)], by="ID_Exoma") 
# MODELO radiomic.0 :25 amostras 
modelo2<- merge(modelo2, dados0[,c(1,3)], by="ID_Exoma")

# Ajuste de um modelo de regressão logística multivariada com todas as variaveis
modelo_logistico <- glm(ClinicalOutcome ~ .  ,family  = binomial(link = 'logit'),  maxit = 1000,data = modelo2[,-c(1)])
summary(modelo_logistico)

# FOREST PLOT ####
print(forest_model(modelo_logistico)+
        theme(text = element_text(size = 14)))

# Encontrando variaveis altamente correlacionadas
vif(modelo_logistico)

#==============================================================================#
# Fazendo a regressao logistica simples do MATH_score classificado
modelo_logistico<-glm(ClinicalOutcome ~ MATH_score , family = binomial(link = 'logit'), data = modelo2[,-c(1)])
summary(modelo_logistico)
# FOREST PLOT ####
print(forest_model(modelo_logistico)+
        theme(text = element_text(size = 14)))

res<-manova(cbind(MATH_score, ADC1, Entropia.Médio) ~ ClinicalOutcome, data=modelo2)
summary(res)




data_matrix <- sampleTableCount.Clinical[,c(2:7,19,24,40)]

# Perform Principal Component Analysis (PCA)
pca_result <- prcomp(data_matrix, center = TRUE, scale = TRUE)

# Extract principal components
pcs <- pca_result$x

# Create a data frame with principal components
pc_data <- data.frame(PC1 = pcs[, 1], PC2 = pcs[, 2])

# Print the summary of the PCA
summary(pca_result)

# Plot the first two principal components
ggplot(pc_data, aes(x = PC1, y = PC2)) +
  geom_point() +
  labs(title = "PCA Plot", x = "Principal Component 1", y = "Principal Component 2")

# Biplot to visualize variable loadings
biplot(pca_result, scale = 0)

# Scree plot to visualize variance explained by each principal component
scree_plot <- ggplot() +
  geom_bar(stat = "identity", aes(x = seq_along(pca_result$sdev), y = pca_result$sdev^2 / sum(pca_result$sdev^2))) +
  labs(title = "Scree Plot", x = "Principal Component", y = "Proportion of Variance Explained")

print(scree_plot)

# Cumulative variance explained plot
cumulative_variance <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)

cumulative_plot <- ggplot() +
  geom_line(aes(x = seq_along(pca_result$sdev), y = cumulative_variance)) +
  labs(title = "Cumulative Variance Explained", x = "Number of Principal Components", y = "Cumulative Variance Explained")

print(cumulative_plot)

concordance_matrix <- matrix(runif(25), ncol = 5)  # Replace this with your actual concordance matrix

# Create a heatmap
heatmap(concordance_matrix, 
        Rowv = NA, Colv = NA,  # Turn off clustering
        col = colorRampPalette(c("white", "blue"))(50),  # Choose a color palette
        main = "Concordance Heatmap",
        xlab = "Samples/Conditions",
        ylab = "Samples/Conditions")

#===============================================================================
# Extract variable loadings
loadings <- as.data.frame(pca_result$rotation[, 1:2])

# Create a custom loading plot using ggplot2
library(ggplot2)

loading_plot <- ggplot(loadings, aes(x = PC1, y = PC2)) +
  geom_point() +
  geom_text(aes(label = rownames(loadings)), hjust = 1, vjust = 1, nudge_x = 0.01) +
  labs(title = "Loading Plot", x = "Principal Component 1", y = "Principal Component 2")

print(loading_plot)

#===============================================================================
