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
#==============================================================================#
# Estamos avaliando os pacientes que RESPONDEM ao tratamento ####
#==============================================================================#
dados0$ClinicalOutcome <-as.factor(ifelse(dados0$Response1 == "nCRT-NR", 0, 1)) 

#==============================================================================#
# Selecionando principais colunas
dados0 <- subset(dados0, 
                select=c("ClinicalOutcome","MATH_score", "Treatment.Protocol","Baseline.T",
                         "CEA.inicial", "Baseline.N","Final.Clinical.Stage","Avoid.TME..1.yes..2.no.",
                         "Distance.from.anal.verge","Gender", "Age", "Tumor.size","MATH_score."))
dados0$ClinicalOutcome<-as.factor(dados0$ClinicalOutcome)
dados0$Tumor.size<-as.numeric(dados0$Tumor.size)
dados0$Baseline.T<-as.factor(dados0$Baseline.T)
dados0$Baseline.T. <-as.factor(ifelse(dados0$Baseline.T == "2", "2", "3-4"))
dados0$Baseline.N<-as.factor(dados0$Baseline.N)
dados0$Final.Clinical.Stage<-as.factor(dados0$Final.Clinical.Stage)
dados0$Final.Clinical.Stage. <-as.factor(ifelse(dados0$Final.Clinical.Stage == "3", "3", "1-2"))
dados0$Avoid.TME..1.yes..2.no.<-as.factor(dados0$Avoid.TME..1.yes..2.no.)
#modelo1$Avoid.TME..1.yes..2.no. <- gsub("2", "0", modelo1$Avoid.TME..1.yes..2.no.)
colunas_caractere <- sapply(dados0, is.character)
# Loop para converter colunas de caractere para fator
for (coluna in names(dados0)[colunas_caractere]) {
  dados0[[coluna]] <- as.factor(dados0[[coluna]])
}
#==============================================================================#
#  DISPERSSÃO DAS MATH_SCORE DAS AMOSTRAS (R/NR) ####
#==============================================================================#
ggplot(dados0, aes(x=MATH_score, y=ClinicalOutcome)) + 
  geom_point() + 
  stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)

#==============================================================================#
# REGRESSÃO LOGISTICA SIMPLES   ####
#==============================================================================#
modelo_logistico<-glm(ClinicalOutcome ~ MATH_score, family= binomial(link = 'logit'), data = dados0)
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
#==============================================================================#
# SELECIONANDO VARIAVEIS COM INFORMAÇOES COMPLETAS PARA AS 52 AMOSTRAS
modelo3<- subset(dados0, 
                  select=c("ClinicalOutcome","MATH_score", 
                           "Treatment.Protocol","Baseline.T","Baseline.N","Gender", "Age",
                           "Distance.from.anal.verge","Tumor.size",
                           "CEA.inicial"))
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



#==============================================================================#

#==============================================================================#
# Fazendo a regressao logistica simples do MATH_score classificado
modelo_logistico<-glm(ClinicalOutcome ~ MATH_score. , family = binomial(link = 'logit'), data = dados0)
summary(modelo_logistico)
# FOREST PLOT ####
print(forest_model(modelo_logistico)+
        theme(text = element_text(size = 14)))


#==============================================================================#

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

