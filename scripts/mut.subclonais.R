
library(data.table)
library(dplyr)
library(reshape2)
library(tidyr)
library(stringr)
library(ggplot2)
library(devtools)
library(dndscv)
library(ggpubr)
library(xlsx)
library(maftools)
library(pheatmap)
library(pROC)
library(tidyverse)

wd <- "D:/PROJETOS-HSL_BP/resultados_Mutect2-2024/"
# wd <- "/media/vande/HDD/PROJETOS-HSL_BP/resultados_Mutect2-2023/"
setwd(wd)
dir.create("subclonais.52samples_subclonal.2024_02_05")
output <- "subclonais.52samples_subclonal.2024_02_05"

input <-  "D:/PROJETOS-HSL_BP/resultados_Mutect2-2023/"


color.R<-"#08415c"
color.NR<-"#cc2936"

#==============================================================================#
#  CARREGANDO TABELA COM TODAS AS MUTAÇÕES
#==============================================================================#
all.mut<-as.data.frame(fread("D:/PROJETOS-HSL_BP/resultados_Mutect2-2024/output-Mutect2_SomaticFilters.2024-02-05/mutation.popFilt.shared_cosmic.tsv", header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))
length(unique(all.mut[,"index"]))  # 56786


all<-as.data.frame(fread("D:/PROJETOS-HSL_BP/resultados_Mutect2-2023/output-Mutect2_SomaticFilters.cov0-maf0.2024-02-05/mutation-ALL.tsv", header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))
#==============================================================================#

#==============================================================================#
#                    STEP 4 - APPLY SOMATIC FILTER        ####
#==============================================================================#

## Quantifica MUTAÇOES not-rare gnomAD
#==============================================================================#
not_rare2exclude <- which(all$freq_gnomAD_lower01=="not-rare")
length(unique(all[not_rare2exclude,"index"]))   # 89

## Quantifica MUTAÇOES not-rare no Abraom
#==============================================================================#
not_rareAbraom2exclude <- which(all$freq_abraom_status=="not-rare")
length(unique(all[not_rareAbraom2exclude,"index"])) #205 

## FILTRO: REMOVE mutações not-rare no gnomAD ou Abraom ####
#==============================================================================#
all.x1 <- all[-c(not_rare2exclude,not_rareAbraom2exclude),]
length(unique(all.x1$index)) 
fwrite(all.x1, paste0(output,"/mutation.popFilt.tsv"), quote = F, sep="\t")
# all.x1<-as.data.frame(fread(paste0(output,"mutation.popFilt.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))

#==============================================================================#
#    EXCLUIR MUTAÇÕES RARAS       ####
#==============================================================================#
## FILTRO: REMOVE mutações raras com <= 3 ocorencias no COSMIC  ####
#==============================================================================#
length(unique(all.x1$index))
rareGnomAD2exclude <- which(all.x1$freq_gnomAD_lower01=="rare" & all.x1$COSMIC_count <= 3)
length(unique(all.x1[rareGnomAD2exclude,"index"])) #excluidas
# View(all.x1[rareGnomAD2exclude,])
rareAbraom2exclude <- which(all.x1$freq_abraom_status=="rare"  & all.x1$COSMIC_count <= 3)
length(unique(all.x1[rareAbraom2exclude,"index"])) #excluidas
# View(all.x1[rareAbraom2exclude,])

all.x2 <- all.x1[-c(rareGnomAD2exclude,rareAbraom2exclude),]
length(unique(all.x2[shared2exclude,"index"])) 



## FILTRO: REMOVE mutações compartilhadas >=2 e  ####
#==============================================================================#
shared2exclude <- which( all.x2$MUTATION_shared >= 2 & all.x2$class != "driver" )
length(unique(all.x2[shared2exclude,"index"]))  
all.xg <- all.x2[-shared2exclude,]  
length(unique(all.xg[,"index"]))  

fwrite(all.xg, paste0(output,"/mutation.popFilt.shared_cosmic.tsv"), quote = F, sep="\t")
all.mut <- data.frame(all.xg)



length(unique(all.mut[,"index"]))  # 56786







#==============================================================================#
# Carregando tabela de dados Clinicos.52 e filtrando mutação dessas amostras ####
#==============================================================================#
Clinical.52 <- read.xlsx("D:/PROJETOS-HSL_BP/Dados Clinicos 2024 ALLversions WXS-TCRseq.xlsx", sheetIndex = "52samples", header = T)
all.mut<- subset(all.mut,SAMPLE %in% Clinical.52$Sample)
length(unique(all.mut[,"index"]))  # 43997

#==============================================================================#
# MARCANDO MUTAÇÕES COM 5+ OCORRENCIAS NO COSMIC
#==============================================================================#
all.mut$COSMIC_count[all.mut$COSMIC_count>= 5] <- paste(all.mut$COSMIC_count[all.mut$COSMIC_count>= 5], "+")

#==============================================================================#
#  QUANTIFICANDO MUTAÇÕES DO MUTECT2 ANTES DE APLICAR FILTROS POPULATIONAL  ####
#==============================================================================#
#  Quantify mutation in HLA genes
#==============================================================================#
genesHLA <- fread("D:/PROJETOS-HSL_BP/resultados_Mutect2-2023/input_Somatic-Filter/genesHLA.txt", col.names = "gene", header = F) 
all.hla<- subset(all.mut, (Gene.refGene %in% genesHLA$gene))
length(unique(all.hla$index))    # 27

# Quantificando regiões das mutações
#==============================================================================#
aux<-(unique(all.mut[,c("index", "Func.refGene")]))  
as.data.frame(table(aux$Func.refGene))

# Quantificando tipos de mutações
#==============================================================================#
aux<-(unique(all.mut[,c("index", "ExonicFunc.refGene")]))  
as.data.frame(table(aux$ExonicFunc.refGene))

# Quantificando Distribuição da MAF 
#==============================================================================#
aux<-(unique(all.mut[,c("index", "MAF")]))  
summary(aux$MAF)
aux1<-as.data.frame(table(all.mut$MAF), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
p<-ggplot(data=aux1, aes(x=Var1, y=Freq))+
  geom_col(width = .005) + 
  geom_vline(xintercept= 0.2, linetype="dashed", color = "red")+ theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
gridExtra::grid.arrange(p, bottom="MAF distribution, from 52 samples")




# Quantificando Distribuição da COV 
#==============================================================================#
aux<-(unique(all.mut[,c("index", "COV")]))  
##PLOT distribuição da COV
#==============================================================================#
aux1<-as.data.frame(table(all.mut$COV), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
p<-ggplot(data=head(aux1, 500), aes(x=Var1, y=Freq))+
  geom_col() + 
  geom_vline(xintercept= 30, linetype="dashed", color = "red")+theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
gridExtra::grid.arrange(p, bottom="COV distribution, from 85 samples")

#==============================================================================#
#                           APLICANDO FILTROS         ####
#==============================================================================#
## FILTRO: REMOVE MUTACOES EM HLA ####
#==============================================================================#
genesHLA <- fread("D:/PROJETOS-HSL_BP/resultados_Mutect2-2023/input_Somatic-Filter/genesHLA.txt", col.names = "gene", header = F) 
all.hla<- subset(all.mut, (Gene.refGene %in% genesHLA$gene))
length(unique(all.hla$index))    # 27

hla2exclude <- which(all.mut$Gene.refGene %in% genesHLA$gene) 
length(unique(all.mut[hla2exclude,"index"]))  #  27
all.mut2 <- all.mut[-hla2exclude,]
length(unique(all.mut2$index))  # 43970


aux<-(unique(all.mut2[,c("index", "MAF")]))
aux1<-as.data.frame(table(all.mut$MAF), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
p<-ggplot(data=aux1, aes(x=Var1, y=Freq))+
  geom_col(width = .005) + 
  geom_vline(xintercept= 0.2, linetype="dashed", color = "red")+ theme_classic()+
  geom_vline(xintercept= 0.15, linetype="dashed", color = "red")+ theme_classic()+
  geom_vline(xintercept= 0.1, linetype="dashed", color = "red")+ theme_classic()+
  geom_vline(xintercept= 0.05, linetype="dashed", color = "red")+ theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
gridExtra::grid.arrange(p, bottom="MAF distribution, from 52 samples")
length(unique( all.mut2[which(all.mut2$MAF<0.2),"index"]))
length(unique( all.mut2[which(all.mut2$MAF<0.15),"index"]))
length(unique( all.mut2[which(all.mut2$MAF<0.1),"index"]))
length(unique( all.mut2[which(all.mut2$MAF<0.05),"index"]))

#==============================================================================#
# Calcular a mediana agrupada por "grupo"
meanMAF <- aggregate(MAF ~ SAMPLE, data = all.mut, FUN = mean)

#==============================================================================#
# EXCLUI AMOSTRAS CUJA MEDIA DAS MAF É ABAIXO DE 0.2
all.mut3<- subset(all.mut2, SAMPLE %in%  meanMAF[meanMAF$MAF>=0.2, ]$SAMPLE)
all.mut.maf<-all.mut3

## FILTRO: REMOVE MUTACOES  MINOR ALLELIC FREQ. < 0.20 ####
#==============================================================================#
maf2exclude <- which(all.mut2$MAF<0.2)
length(unique(all.mut2[maf2exclude,"index"]))  #  27136  20%:20522; 15%:15271; 10%:11445
all.mut.maf <- all.mut2[-maf2exclude,]
length(unique(all.mut.maf$index))  #29659



## SELECIONANDO AS MUTAÇÕES QUE DE BAIXO MAF ####
#==============================================================================#
all.mut.maf.excluded <- all.mut2[maf2exclude,]
length(unique(all.mut.maf.excluded$index))  # 20522
aux1<-as.data.frame(table(all.mut.maf$MAF), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
p<-ggplot(data=aux1, aes(x=Var1, y=Freq))+
  geom_col(width = .005) + 
  geom_vline(xintercept= 0.05, linetype="dashed", color = "red")+ theme_classic()+
  geom_vline(xintercept= 0.1, linetype="dashed", color = "red")+ theme_classic()+
  geom_vline(xintercept= 0.15, linetype="dashed", color = "red")+ theme_classic()+
  geom_vline(xintercept= 0.2, linetype="dashed", color = "red")+ theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
gridExtra::grid.arrange(p, bottom="MAF distribution, from 52 samples")

# all.mut.maf<-all.mut.maf.excluded

# ## FILTRO: REMOVE MUTACOES COV < 30 ####
# #==============================================================================#
cov2exclude <- which(all.mut.maf$COV<30)
length(unique(all.mut.maf[cov2exclude,"index"]))  #  431
all.mut.cov <- all.mut.maf[-cov2exclude,]
length(unique(all.mut.cov$index))  # 23044

## FILTRO: SELECT EXONIC AND SPLICING ####
#==============================================================================#
table(all.mut.cov$Func.refGene)
all.Coding <- subset(all.mut.cov, Func.refGene %in% c("splicing", "exonic", "exonic;splicing", "splicing;exonic"))
length(unique(all.Coding$index))  # 22937
unique(all.Coding$Func.refGene)
#34126765 - xGen Exome Research Panel V2 Target Regions
Mb <- 34126765/1000000    
#Mb=34.126
tmb <- as.data.frame(table(all.Coding$SAMPLE))
names(tmb) <- c("SAMPLE", "N_Mutations_TMB")
tmb$TMB <- round(tmb$N_Mutations/Mb,2)

fwrite(tmb, paste0(output,"/TMB.tsv"), quote = F, sep="\t")
fwrite(all.Coding, paste0(output,"/mutation_codRegion_Freq005.tsv"), quote = F, sep="\t")

#==============================================================================#
#all.Coding<-as.data.frame(fread(paste0(output,"mutation_codRegion_Freq005.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))
aux<-(unique(all.Coding[,c("index", "ExonicFunc.refGene")]))  
as.data.frame(table(aux$ExonicFunc.refGene))

#  REMOVE MUTACOES CLASSIFICADAS COM SINONIMAS (synonymous SNV)   11218
#==============================================================================#
somatic = subset(all.Coding, !(ExonicFunc.refGene %in% c("synonymous SNV")))
length(unique(somatic[,"index"]))     # 17167





### SALVAR A TABELA COM MUTACOES SOMATICAS  ####
#==============================================================================#
fwrite(somatic, paste0(output,"/somatic_mutation.tsv"), quote = F, sep="\t")
#==============================================================================#
#somatic<-as.data.frame(fread(paste0(output,"/somatic_mutation.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))

#==============================================================================#
# CALCULAR MATH SCORE  ####
#==============================================================================#
somatic.LENGTH <- as.data.frame(tapply(somatic$MAF, somatic$SAMPLE, length))
somatic.LENGTH$SAMPLE <- row.names(somatic.LENGTH)
names(somatic.LENGTH) <- c("N_mutations", "SAMPLE")
somatic.MEDIAN <- as.data.frame(tapply(somatic$MAF, somatic$SAMPLE, median))
somatic.MEDIAN$SAMPLE <-row.names(somatic.MEDIAN)
names(somatic.MEDIAN) <- c("MAF_median", "SAMPLE")
somatic.MAD <- as.data.frame(tapply(somatic$MAF, somatic$SAMPLE, mad))
somatic.MAD$SAMPLE <-row.names(somatic.MAD)
names(somatic.MAD) <- c("MAF_MAD", "SAMPLE")
match_somatic <- match(somatic.MAD$SAMPLE, somatic.MEDIAN$SAMPLE)
somatic.MAD$MAF_median <- somatic.MEDIAN$MAF_median[match_somatic]
match_somatic_length <- match(somatic.MAD$SAMPLE, somatic.LENGTH$SAMPLE)
somatic.MAD$N_mutations <- somatic.LENGTH$N_mutations[match_somatic_length]
somatic.MATH <- somatic.MAD
somatic.MATH$MATH_score <- 100*(somatic.MATH$MAF_MAD/somatic.MATH$MAF_median)

#==============================================================================#
# TABELA COM SCORES DE HETEROGENIEDADE  ####
#==============================================================================#
somatic_Samples<- merge(somatic.MATH, tmb, by = "SAMPLE")

# SALVAR A TABELA MATH-SCORE, TMB DE CADA AMOSTRA  ###
#==============================================================================#
fwrite(somatic_Samples, paste0(output,"/MATH_TMB-SAMPLES.tsv"), quote = F, sep="\t")

#==============================================================================#
# Carregando tabela de dados Clinicos ####
#==============================================================================#
Clinical.52 <- read.xlsx("D:/PROJETOS-HSL_BP/Dados Clinicos 2024 ALLversions WXS-TCRseq.xlsx", sheetIndex = "52samples", header = T)
Clinical.52<- Clinical.52[rowSums(is.na(Clinical.52)) != ncol(Clinical.52), -c(6)]


# SALVAR ARQUIVOS BED PARA FILTRAR VCF ####
#==============================================================================#
dir.create(paste0(output,"mutations_to_filter_fromVCF/"))
samples <-unique(Clinical.52$Sample)
for (sp in samples) {
  BED<-subset(somatic, SAMPLE==sp, select = c("CHR", "POS","REF","ALT","index"))
  BED<-unique(BED)
  fwrite(BED,paste0(output,"mutations_to_filter_fromVCF/",sp,".bed"), sep = "\t", col.names = T, row.names = F, quote = F)
}
# SALVAR ARQUIVOS BED COM SINONIMAS PARA FILTRAR VCF ####
#==============================================================================#
samples<-unique(Clinical.52$Sample)
dir.create(paste0(output,"mutations_to_filter_fromVCF/"))
for (sp in samples) {
  BED<-subset(all.Coding, SAMPLE==sp, select = c("CHR", "POS","REF","ALT","index"))
  BED<-unique(BED)
  fwrite(BED,paste0(output,"mutations_to_filter_fromVCF/",sp,".synonymous.bed"), sep = "\t", col.names = T, row.names = F, quote = F)
}

## Carregando dados de MATH_SCORE, TMB, MUT, SAMPLE
sampleTableCount<- fread(paste0(output,"/MATH_TMB-SAMPLES.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)


#==============================================================================#
## JUNTANDO TABELAS COUNT_MUT_MATH_TMB +  dados_clinicos 
#==============================================================================#
sampleTableCount.Clinical <- merge(sampleTableCount, Clinical.52, by.x = "SAMPLE", by.y = "Sample")
sampleTableCount.Clinical2 <- subset(sampleTableCount.Clinical, !(ID_Exoma %in% c("ROP-107", "ROP-83")), )

#==============================================================================#
### PLOT MATH_score ####
#==============================================================================#
## MATH-SCORE colour by ClinicalOutcome 

# PLOT BOXPLOT MATH_score CORRELATION CLINICAL
#==============================================================================#
# ClinicalOutcome
table(sampleTableCount.Clinical$ClinicalOutcome)
aux<-(table(sampleTableCount.Clinical$ClinicalOutcome))
p<-ggplot(data=sampleTableCount.Clinical, aes(x= reorder(ID_Exoma, -MATH_score), y=MATH_score, fill=ClinicalOutcome)) +
  geom_bar(stat="identity",color="black")+theme_classic() +
  scale_fill_manual(values = c(color.R, color.NR), breaks = c("nCRT-R", "nCRT-NR")) +
  geom_hline(yintercept= median(sampleTableCount.Clinical$MATH_score), linetype="dashed", color = "black")+
  annotate(geom="text", x=30, y=50, label=paste("median MATH_score= ",sprintf("%.2f",median(sampleTableCount.Clinical$MATH_score))), color="black")+
  xlab("SAMPLES") +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
gridExtra::grid.arrange(p, bottom=paste0("nCRT-R= ",aux[2], "; nCRT-NR= ",aux[1]))

shapiro.test(sampleTableCount.Clinical$MATH_score)
p <- ggplot(sampleTableCount.Clinical  , aes(x=ClinicalOutcome, y=MATH_score, fill=ClinicalOutcome)) + 
  geom_boxplot()+ 
  scale_fill_manual(values = c(color.R, color.NR), breaks = c("nCRT-R", "nCRT-NR")) +
  annotate("text", label = paste('N = ',aux[2], sep = ''), x = 2, y = 3, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = ',aux[1], sep = ''), x = 1, y = 3, angle = 0, size = 3.75) +
  theme_classic()+theme(legend.position="none")
table(sampleTableCount.Clinical$ClinicalOutcome)
result_test <- t.test(MATH_score ~ ClinicalOutcome, data = sampleTableCount.Clinical)
gridExtra::grid.arrange(p, bottom= paste(result_test$method,"; p-value: ",sprintf("%.4f",result_test$p.value)))
result_test


#==============================================================================#
#                   DISTRIBUICAO DAS MAF NAS 52 AMOSTRAS  ####
#==============================================================================#
# Para plotar a distuibuiC'C#o dos MAFs estou usando as mutaC'C5es NON_SYNONIMOUS + SINONIMOUS
#somatic.52<-as.data.frame(fread(paste0(output,"/somatic_mutation.72samples.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))

aux<-as.data.frame(table(all.mut.maf.excluded$SAMPLE))
summary(aux$Freq)
hist(aux$Freq)
boxplot(aux$Freq)

auxx<-somatic
aux1 <- dplyr::select(auxx, SAMPLE, MAF)
aux2 <- as.data.frame(table(aux1), stringsAsFactors = F)
aux2$MAF <- as.double(aux2$MAF)
aux2 <- merge(aux2, sampleTableCount.Clinical[,c("SAMPLE","ID_Exoma","ClinicalOutcome")], by= "SAMPLE", all.x = T)

# MAFs das amostras outliers    ####
#==============================================================================#
# # nCRT-R ClinicalOutcome 
P<- ggplot(data=aux2,
           aes(x=MAF, y=Freq, group = ID_Exoma, colour= ClinicalOutcome))+
  geom_col() +
  scale_y_sqrt()+
  # scale_x_continuous(breaks=seq(0,1,0.1))+
  facet_wrap(vars(ID_Exoma))+theme_minimal()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")
gridExtra::grid.arrange(P, bottom="MAF distribution from Excluded Mutations (52 samples)")

