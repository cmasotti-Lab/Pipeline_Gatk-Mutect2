library("forestmodel")
library(car)

# CARREGANDO DADOS DA TABELA dados GENADOS NO SCRIPT ANTERIOR
# Organizando os dados continuos em (numeric) 
# E os dados categoricos em (factor)
dados0<-dados

#dados0$ClinicalOutcome <-ifelse(dados0$ClinicalOutcome == "nCRT-R", 1, 0) # USAR VITAL.STATUS OU RESPOSTA, RESPOSTA2
dados0$CEA.inicial<-as.numeric(dados0$CEA.inicial)
dados0$Clones<-as.numeric(dados0$Clones)
dados0$Tumor.size<-as.numeric(dados0$Tumor.size)
dados0$Baseline.T<-as.factor(dados0$Baseline.T)
dados0$Baseline.N<-as.factor(dados0$Baseline.N)
dados0$Final.Clinical.Stage<-as.factor(dados0$Final.Clinical.Stage)
dados0$Avoid.TME..1.yes..2.no.<-as.factor(dados0$Avoid.TME..1.yes..2.no.)
#cenario1$Avoid.TME..1.yes..2.no. <- gsub("2", "0", cenario1$Avoid.TME..1.yes..2.no.)
# Identificar colunas de caractere e transformando em factor
colunas_caractere <- sapply(dados0, is.character)
# Loop para converter colunas de caractere para fator
for (coluna in names(dados0)[colunas_caractere]) {
  dados0[[coluna]] <- as.factor(dados0[[coluna]])
}
colnames(dados0)[colnames(dados0) == "Tratamento"] <- "Treatment.Protocol"

#==============================================================================#
# Fazendo a regressao liner bivariadas, ou seja, uma variante por vez
modelo_bivar<- summary(lm(ClinicalOutcome ~ MATH_score, data = dados0))
modelo_bivar
print(forest_model(modelo_bivar)+
  theme(text = element_text(size = 14)))
#==============================================================================#
# Selecionando principais colunas
dados0<- subset(dados0, 
                select=c("ClinicalOutcome","MATH_score", "Treatment.Protocol","Baseline.T",
                         "CEA.inicial", "Baseline.N","Final.Clinical.Stage","Avoid.TME..1.yes..2.no.",
                         "Distance.from.anal.verge","Gender", "Age", "Tumor.size","MATH_score."))

#==============================================================================#
# CENARIO 1: 52 AMOSTRAS
cenario1<- subset(dados0, 
                select=c("ClinicalOutcome","MATH_score", 
                         "Treatment.Protocol","Baseline.T","Baseline.N","Gender", "Age"))

# Ajuste de um modelo de regressão logística multivariada com todas as variaveis
modelo_logistico <- glm(ClinicalOutcome ~ .  ,family  = binomial(link = 'logit'), data = cenario1)
summary(modelo_logistico)

# FOREST PLOT ####
print(forest_model(modelo_logistico)+
        theme(text = element_text(size = 14)))

# Encontrando variaveis altamente correlacionadas
 vif(modelo_logistico)
#==============================================================================#
#==============================================================================#
# CENARIO 2
cenario2<- subset(dados0, 
                  select=c("ClinicalOutcome","MATH_score", 
                           "Treatment.Protocol","Baseline.T","Baseline.N","Gender", "Age",
                           "Final.Clinical.Stage","Distance.from.anal.verge","Tumor.size"))

# Ajuste de um modelo de regressão logística multivariada com todas as variaveis
modelo_logistico <- glm(ClinicalOutcome ~ . ,  data = cenario2)
summary(modelo_logistico)
vif(modelo_logistico)
# CENARIO 2
var.cenario2<-c("ClinicalOutcome","MATH_score","Treatment.Protocol","Baseline.T",
                "Baseline.N","Gender", "Age","Distance.from.anal.verge","Tumor.size")

cenario2<- subset(dados0, 
                  select=c("ClinicalOutcome","MATH_score", 
                           "Treatment.Protocol","Baseline.T","Baseline.N","Gender", "Age",
                           "Distance.from.anal.verge","Tumor.size"))

# Ajuste de um modelo de regressão logística multivariada com todas as variaveis
modelo_logistico <- glm(ClinicalOutcome ~ . , data = cenario2)
summary(modelo_logistico)

# FOREST PLOT ####
print(forest_model(modelo_logistico)+   
        theme(text = element_text(size = 14)))

# Encontrando variaveis altamente correlacionadas
 vif(modelo_logistico)
#==============================================================================#
#==============================================================================#
# SELECIONANDO VARIAVEIS COM INFORMAÇOES COMPLETAS PARA AS 52 AMOSTRAS
cenario3<- subset(dados0, 
                  select=c("ClinicalOutcome","MATH_score", 
                           "Treatment.Protocol","Baseline.T","Baseline.N","Gender", "Age",
                           "Distance.from.anal.verge","Tumor.size",
                           "CEA.inicial"))

# Ajuste de um modelo de regressão logística multivariada com todas as variaveis
modelo_logistico <- lm(ClinicalOutcome ~ .  , data = cenario3)
summary(modelo_logistico)

# FOREST PLOT ####
print(forest_model(modelo_logistico)+
        theme(text = element_text(size = 14)))

# Encontrando variaveis altamente correlacionadas
vif(modelo_logistico)
#==============================================================================#



var.cenario2 <- c("MATH_score", "Treatment.Protocol", "Baseline.T",
                  "Baseline.N", "Gender", "Age", "Distance.from.anal.verge", "Tumor.size",
                  "Final.Clinical.Stage", "CEA.inicial")

results <- data.frame(Variable = character(), P_value = numeric(), stringsAsFactors = FALSE)

for (variable in var.cenario2) {
  formula <- as.formula(paste("ClinicalOutcome ~", variable))
  modelo_bivar <- lm(formula, data = dados0)
  
  # Extract the p-value from the summary
  p_value <- summary(modelo_bivar)$coefficients[, "Pr(>|t|)"][2]
  print(summary(modelo_bivar))
  # Store results in a data frame
  results <- rbind(results, data.frame(Variable = variable, P_value = p_value))
}

# Print the formatted table
knitr::kable(results)





#==============================================================================#
#==============================================================================#
# Fazendo a regressao liner bivariadas, ou seja, uma variante por vez
modelo_logistico<-summary(lm(ClinicalOutcome ~ MATH_score., data = dados0))
modelo_logistico
# FOREST PLOT ####
print(forest_model(modelo_logistico)+
        theme(text = element_text(size = 14)))












