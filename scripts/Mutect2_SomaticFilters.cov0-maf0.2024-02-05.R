### VLIRA FEV2024 ####
### Pipeline de indentificacao de mutacoes somaticas usando Mutect2 ###
# ==============================================================================#
# Na versão anterior (2023-09-26) faltava excluir mutações na lista de HLA
# Apliquei os filtros:
# gnomAD: 0.005
# Abraom: 0.005
# shared >=2 & cosmic <=3
# MAF < 0.20
# COV < 30
# Region: exonic, splincing
# Excluido outliers: ROP-129, ROP-91, ROP-81
# Em alguns momentos é removido: ROP-107 e ROP-83
# Inserido legendas nas Figuras
# Refinando Calculo do dNdS 
#   Fiz 4 teste com as colunas (Chr,Start,Ref,Alt) ao invéz de (CHR,POS,REF,ALT), 
#   e com o genoma referencia GencodeV18 e GRCh38.p12
# MAS AINDA NÃO CONSEGUI CORRIGIR O DNDS_CV PARA OLHAR APENAS AS TARGETS DO EXOMA
# Decidi usar o dndsout2,  (Chr,Start,Ref,Alt) e GencodeV18
# ONCOPLOT agora mostras as mutaçoes sinonimas

# Fisher Test (dnds_genes, dMMR, drivers)
# Calcula Curva ROC do TMB e MATH_score para Resposta1

#==============================================================================#

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

wd <- "D:/PROJETOS-HSL_BP/resultados_Mutect2-2023/"
# wd <- "/media/vande/HDD/PROJETOS-HSL_BP/resultados_Mutect2-2023/"
setwd(wd)
dir.create("output-Mutect2_SomaticFilters.cov0-maf0.2024-02-05")
output <- "output-Mutect2_SomaticFilters.cov0-maf0.2024-02-05/"

#==============================================================================#
# STEP 0 -  Datasets reorganization ####
#==============================================================================#
snps<- fread(file="input_Somatic-Filter/annovar.norm.hg38_multianno.txt", header=T, sep="\t", na.strings=c(".","NA"), fill=TRUE, check.names = FALSE)
snps<- as.data.frame(snps)
samples_ID <-fread("input_Somatic-Filter/sampleID.list", header = F)

# Selecting main columns
#==============================================================================#
snps.flt<-  dplyr::select(snps, c(1:5, ends_with('.refGene'),"gnomAD_exome_ALL",starts_with("abraom"),"cosmic95", starts_with("ICGC"),"Interpro_domain", starts_with("Otherinf") ))
colnames(snps.flt) <-  c( colnames(snps.flt[,c(1:23)]),"Otherinfo1","Otherinfo2","Otherinfo3","CHR","POS","ID","REF","ALT","QUAL","FILTER", "INFO", "FORMAT", 
                          samples_ID$V1)

# Creating INDEX (CHR,POS,REF,ALT) and select the main columns
#==============================================================================#
snps.flt$index <- paste(snps.flt$CHR, snps.flt$POS, snps.flt$REF, snps.flt$ALT, sep="_")
exCol=colnames(snps.flt[,-which(names(snps.flt) %in% c("Otherinfo1","Otherinfo2","Otherinfo3","ID","index"))])
snps.flt <-  dplyr::select(snps.flt, index,all_of(exCol))
dim(snps.flt)   #371253 
fwrite(snps.flt, paste0(output,"/tab0.tsv"), quote = F, sep="\t")

#==============================================================================#
#                 STEP 1 - Clear dataset ####
#==============================================================================#
# tab1<-as.data.frame(fread(paste0(output,"/tab0.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))
tab1<- snps.flt

## Select mutations "PASS" from Mutect2 ###
#==============================================================================#
# salva as NOT PASS em uma tabela
fwrite(subset(tab1, FILTER != "PASS"), paste0(output,"/mutation-ALL-notPASS.tsv"), quote = F, sep="\t")
tab1 = subset(tab1, FILTER == "PASS")
dim(tab1)   #59957            

## Restructure tables, all samples in the same column ###
#==============================================================================#
tab1 <- reshape2::melt(tab1,id=colnames(tab1[,1:which(names(tab1)== "FORMAT")]))

## Remove homozygous alleles to REF ###
#==============================================================================#
tab1 <- tab1[!grepl("^./.:.:.:", tab1$value),]
tab1 <- tab1[!grepl("^0/0", tab1$value),]
tab1 <- tab1[!grepl("^0\\|0", tab1$value),]
fwrite(tab1, paste0(output,"/tab1.tsv"), quote = F, sep="\t")

#==============================================================================#
#             STEP 2 - Restructure and filter the tables ####
#==============================================================================#
# tab1<-as.data.frame(fread(paste0(output,"/tab1.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))
length(unique(tab1$index))   #59957         

## Separate genotype info ("GT","AD","AF","DP","F1R2","F2R1","FAD") ###
#==============================================================================#
tab2<- separate(data = tab1, col = value, into = c("GT","AD","AF","DP","F1R2","F2R1","FAD"), sep = ":")

## Separate Alleles ###
#==============================================================================#
tab2$aux = tab2$GT
tab2$aux = str_replace(tab2$aux,"/","|")   
tab2<-cbind(tab2, read.table(text = as.character(tab2$aux), sep = "|"))
names(tab2)[c(ncol(tab2)-1,ncol(tab2))] <- c('GT1','GT2')
tab2<- subset(tab2, select= -c(aux))
setDT(tab2)[,paste0("AD", 1:2) := tstrsplit(AD, ",")]

## Clear memory ###
#==============================================================================#
remove(snps, tab1, snps.flt)
gc()

## CALCULAR MAF e COV ###
#==============================================================================#
tab2$COV <- as.numeric(tab2$AD1)+as.numeric(tab2$AD2)
tab2$MAF <- round(as.numeric(tab2$AD2)/(as.numeric(tab2$AD1)+as.numeric(tab2$AD2)),5)
colnames(tab2)[colnames(tab2) == "variable"] <- "SAMPLE"

# Quantificando mutações compartilhadas >=2
#==============================================================================#
aux <- data.frame(table(tab2$index))
tab2 <- merge(tab2, aux, by.x = "index", by.y = "Var1")
colnames(tab2)[names(tab2)=="Freq"]<-"MUTATION_shared"
nrow(unique(tab2[,"index"]))                   #59957
fwrite(tab2, paste0(output,"/tab2.mutation-ALL.tsv"), quote = F, sep="\t")



#==============================================================================#
##  QUANTIFICANDO MUTAÇÕES DO ALL MUTECT ANTES DE APLICAR FILTROS POPULATIONAL  ####
#==============================================================================#
all<-as.data.frame(fread(paste0(output,"tab2.mutation-ALL.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))
length(unique(all[,"index"]))                   #59957

# Quantify mutation in HLA genes
#==============================================================================#
genesHLA <- fread("input_Somatic-Filter/genesHLA.txt", col.names = "gene", header = F) 
all.hla<- subset(all, (Gene.refGene %in% genesHLA$gene))
length(unique(all.hla$index))    #37 mut em HLA

# Quantificando regiões das mutações
#==============================================================================#
aux<-(unique(all[,c("index", "Func.refGene")]))  
as.data.frame(table(aux$Func.refGene))

# Quantificando tipos de mutações
#==============================================================================#
aux<-(unique(all[,c("index", "ExonicFunc.refGene")]))  
as.data.frame(table(aux$ExonicFunc.refGene))

# Quantificando mutações compartilhadas >=2
#==============================================================================#
shared2exclude <- which((all$MUTATION_shared >= 2) )
length(unique(all[shared2exclude,"index"]))  #3048 

# Quantificando mutações compartilhadas >=2, Quantas são em genes drivers? 
#==============================================================================#
dgenes<-fread("input_Somatic-Filter/driver-genes.txt", header = F, col.names = "gene")
shared2exclude <- which((all$MUTATION_shared >= 2) & (all$Gene.refGene %in% dgenes$gene))
length(unique(all[shared2exclude,"index"]))  #44 

# Quantificando Distribuição da COV
#==============================================================================#
aux<-(unique(all[,c("index", "COV")]))  
summary(aux$COV)

# Quantificando Distribuição  MAF 
#==============================================================================#
aux<-(unique(all[,c("index", "MAF")]))  
summary(aux$MAF)

#==============================================================================#
#             STEP 3 - MARK Rare mutations ####
#==============================================================================#
# all<-fread(paste0(output,"tab2.mutation-ALL.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)

## Marking rare mutations from gnomAD and Abraom 
#==============================================================================#
all$freq_gnomAD_lower01[is.na(all$gnomAD_exome_ALL)] <- "not-reported"
all$freq_gnomAD_lower01[all$gnomAD_exome_ALL<0.005] <- "rare"
all$freq_gnomAD_lower01[all$gnomAD_exome_ALL>=0.005] <- "not-rare"
all$freq_abraom_status[is.na(all$abraom_freq)] <- "not-reported"
all$freq_abraom_status[all$abraom_freq<0.005] <- "rare"    
all$freq_abraom_status[all$abraom_freq>=0.005] <- "not-rare"  

## Quantify mutations present on COSMIC 
#==============================================================================#
all <- separate(data = all, col = cosmic95, into = c("COSMIC_ID","OCCURENCE"), sep = "OCCURENCE=")
all$COSMIC_OCCURENCE <- all$OCCURENCE
all <- separate(data =all, col =OCCURENCE, into = c("V_1","V_2","V_3","V_4","V_5"), sep = ",")
aux <- as.data.frame(lapply(all[,c("V_1","V_2","V_3","V_4","V_5")], function(y) as.numeric(gsub("\\(.*\\)", "", y))))
aux[is.na(aux)] <- 0
all <- cbind(subset(all, select = -c(V_1,V_2,V_3,V_4,V_5)), aux)
all$COSMIC_count <- rowSums(all[,c("V_1","V_2","V_3","V_4","V_5")])
all <- select(all, !starts_with("V_"))

## Classify mutations like (none, driver, dMMMR) 
#==============================================================================#
dgenes<-fread("input_Somatic-Filter/driver-genes.txt", header = F, col.names = "gene")
dMMR<-fread("input_Somatic-Filter/dMMR_genes.txt", header = F, col.names = "gene")
all$class = "none"
all$class[all$Gene.refGene %in% dgenes$gene] <-"driver"
all$class[all$Gene.refGene %in% dMMR$gene] <-"dMMR"

#==============================================================================#
# ULTIMA TABALA DE MUTAÇÕES ANTES DE APLICAR FILTRO SOMATICO
#==============================================================================#
fwrite(all, paste0(output,"/mutation-ALL.tsv"), quote = F, sep="\t")
all<-as.data.frame(fread(paste0(output,"mutation-ALL.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))
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

## FILTRO: REMOVE mutações compartilhadas >=2 e <= 3 ocorencias no COSMIC ####
#==============================================================================#
shared2exclude <- which( all.x1$MUTATION_shared >= 2 & all.x1$COSMIC_count <= 3)
length(unique(all.x1[shared2exclude,"index"]))  
all.xg <- all.x1[-shared2exclude,]  
length(unique(all.xg[,"index"]))             
fwrite(all.xg, paste0(output,"/mutation.popFilt.shared_cosmic.tsv"), quote = F, sep="\t")
all.mut <- data.frame(all.xg)






#==============================================================================#
#  CARREGANDO TABELA COM TODAS AS MUTAÇÕES
all.mut<-as.data.frame(fread(paste0(output,"mutation.popFilt.shared_cosmic.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))
#==============================================================================#
#SELECIONANTO APENAS AS 52 AMOSTRAS NAO METASTATICAS
samples.52 <- sampleTableCount.Clinical$SAMPLE
all.mut<- subset(all.mut,SAMPLE %in% samples.52)
length(unique(all.mut[,"index"]))  
#==============================================================================#
all.mut$COSMIC_count[all.mut$COSMIC_count>= 5] <- paste(all.mut$COSMIC_count[all.mut$COSMIC_count>= 5], "+")

#==============================================================================#
#  QUANTIFICANDO MUTAÇÕES DO ALL MUTECT ANTES DE APLICAR FILTROS POPULATIONAL  ####
#==============================================================================#

#  Quantify mutation in HLA genes
#==============================================================================#
genesHLA <- fread("input_Somatic-Filter/genesHLA.txt", col.names = "gene", header = F) 
all.hla<- subset(all.mut, (Gene.refGene %in% genesHLA$gene))
length(unique(all.hla$index))    # 35

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
summary(aux$COV)

##PLOT distribuição da COV
#==============================================================================#
aux1<-as.data.frame(table(all.mut$COV), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
p<-ggplot(data=head(aux1, 500), aes(x=Var1, y=Freq))+
  geom_col() + 
  geom_vline(xintercept= 30, linetype="dashed", color = "red")+theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
gridExtra::grid.arrange(p, bottom="COV distribution, from 85 samples")


## FILTRO: REMOVE MUTACOES EM HLA ####
#==============================================================================#
genesHLA <- fread("input_Somatic-Filter/genesHLA.txt", col.names = "gene", header = F) 
all.hla<- subset(all.mut, (Gene.refGene %in% genesHLA$gene))
length(unique(all.hla$index))    # 35

hla2exclude <- which(all.mut$Gene.refGene %in% genesHLA$gene) 
length(unique(all.mut[hla2exclude,"index"]))  #  35
all.mut2 <- all.mut[-hla2exclude,]
length(unique(all.mut2$index))  #56751

## FILTRO: REMOVE MUTACOES  MINOR ALLELIC FREQ. < 0.20 ####
#==============================================================================#
maf2exclude <- which(all.mut2$MAF<0.2)
length(unique(all.mut2[maf2exclude,"index"]))  #  27136  20%:20522; 15%:15271; 10%:11445
all.mut.maf <- all.mut2[-maf2exclude,]
length(unique(all.mut.maf$index))  #29659

all.mut.maf.excluded <- all.mut2[maf2exclude,]
aux1<-as.data.frame(table(all.mut.maf.excluded$MAF), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
p<-ggplot(data=aux1, aes(x=Var1, y=Freq))+
  geom_col(width = .005) + 
  geom_vline(xintercept= 0.05, linetype="dashed", color = "red")+ theme_classic()+
  geom_vline(xintercept= 0.1, linetype="dashed", color = "red")+ theme_classic()+
  geom_vline(xintercept= 0.15, linetype="dashed", color = "red")+ theme_classic()+
  geom_vline(xintercept= 0.2, linetype="dashed", color = "red")+ theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
gridExtra::grid.arrange(p, bottom="MAF distribution, from 85 samples")

# ## FILTRO: REMOVE MUTACOES COV < 30 ####
# #==============================================================================#
cov2exclude <- which(all.mut.maf$COV<30)
length(unique(all.mut.maf[cov2exclude,"index"]))  #  541
all.mut.cov <- all.mut.maf[-cov2exclude,]
length(unique(all.mut.cov$index))  # 29119

all.mut.cov<- all.mut.maf
length(unique(all.mut.cov$index))

## FILTRO: SELECT EXONIC AND SPLICING ####
#==============================================================================#
all.Coding <- subset(all.mut.cov, Func.refGene %in% c("splicing", "exonic", "exonic;splicing", "splicing;exonic"))
length(unique(all.Coding$index))  #28978
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
length(unique(somatic[,"index"]))     # 21399

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
Clinical <- read.xlsx("Dados Clinicos JP final  - Dez2022 - FILTRADO.xlsx", sheetIndex = 1, header = T)
Clinical<- Clinical[rowSums(is.na(Clinical)) != ncol(Clinical), -c(6)]

## Renomeando desfeichos clinicos ####
#==============================================================================#
Clinical$Response1 <- gsub("Completa", "nCRT-R" , Clinical$Resposta)
Clinical$Response1 <- gsub("Incompleta", "nCRT-NR" , Clinical$Response1)
Clinical$Response2 <- gsub("Completa", "nCRT-R" , Clinical$Resposta2)
Clinical$Response2 <- gsub("Incompleta Sem Metastase", "nCRT-NR Not-Metastatic" , Clinical$Response2)
Clinical$Response2 <- gsub("Incompleta Com Metastase", "nCRT-NR Metastatic" , Clinical$Response2)
Clinical$Response3 <- gsub("nao metastatico", "Not-Metastatic" , Clinical$Resposta3)
Clinical$Response3 <- gsub("metastatico", "Metastatic" , Clinical$Response3)
#==============================================================================#

## Listras de amostras excluidas pelo corpo clinico 
samples.excl <- c("ROP-116","ROP-118","ROP-120","ROP-121","ROP-122","ROP-123","ROP-125","ROP-126",
                  "ROP-127","ROP-130","ROP-131","ROP-132","ROP-23","ROP-84", "ROP-129", "ROP-81", "ROP-91")

## Selecionando amostas do Exoma 
Clinical <- Clinical[Clinical$WXS== "Y",]
## Excluindo amostras excluidas pelo corpo clinico 
Clinical <- subset(Clinical, !(ID_Exoma %in% samples.excl) )


# SALVAR ARQUIVOS BED PARA FILTRAR VCF ####
#==============================================================================#
samples<-unique(Clinical$Sample)
dir.create(paste0(output,"mutations_to_filter_fromVCF/"))
for (sp in samples) {
  BED<-subset(somatic, SAMPLE==sp, select = c("CHR", "POS","REF","ALT","index"))
  BED<-unique(BED)
  fwrite(BED,paste0(output,"mutations_to_filter_fromVCF/",sp,".bed"), sep = "\t", col.names = T, row.names = F, quote = F)
}

# SALVAR ARQUIVOS BED COM SINONIMAS PARA FILTRAR VCF ####
#==============================================================================#
samples<-unique(Clinical$Sample)
dir.create(paste0(output,"mutations_to_filter_fromVCF/"))
for (sp in samples) {
  BED<-subset(all.Coding, SAMPLE==sp, select = c("CHR", "POS","REF","ALT","index"))
  BED<-unique(BED)
  fwrite(BED,paste0(output,"mutations_to_filter_fromVCF/",sp,".synonymous.bed"), sep = "\t", col.names = T, row.names = F, quote = F)
}

## Carregando dados de MATH_SCORE, TMB, MUT, SAMPLE
sampleTableCount<- fread(paste0(output,"/MATH_TMB-SAMPLES.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)

aux<-subset(somatic_Samples, SAMPLE %in% Clinical$Sample)
fwrite(aux, paste0(output,"/MATH_TMB-69.SAMPLES.tsv"), quote = F, sep="\t")

## Filtrando mutações dos 69 pacientes finais 
sampleTableCount<-subset(sampleTableCount, (SAMPLE %in% Clinical$Sample) )

## Juntando tabelas TMB_MATH_#MUT +  dados_clinicos 
sampleTableCount.Clinical <- merge(sampleTableCount, Clinical, by.x = "SAMPLE", by.y = "Sample")


sampleTableCount.Clinical2 <- subset(sampleTableCount.Clinical, !(ID_Exoma %in% c("ROP-107", "ROP-83")), )

## Inserindo coluna ID_Exoma na tabela somatic 
somatic <- merge(somatic, Clinical[,c("Sample","ID_Exoma")], by.x = "SAMPLE", by.y = "Sample", all.x = T)

## Removendo amostras excluidas por criterios clinicos da tabela somatic 
somatic.69 <- somatic[!is.na(somatic$ID_Exoma),] 

### SALVAR A TABELA COM MUTACOES SOMATICAS DAS 69 amostras ####
#==============================================================================#
fwrite(somatic.69, paste0(output,"/somatic_mutation.69samples.tsv"), quote = F, sep="\t")

#==============================================================================#
### PLOT MATH_score ####
#==============================================================================#
sampleTableCount.Clinical<-subset(sampleTableCount.Clinical, Response3 == "Not-Metastatic")
sampleTableCount.Clinical<-subset(sampleTableCount.Clinical, ID_Exoma  != "ROP-50")

## MATH-SCORE colour by Resposta1 
p<-ggplot(data=sampleTableCount.Clinical, aes(x= reorder(ID_Exoma, -MATH_score), y=MATH_score, fill=Response1)) +
  geom_bar(stat="identity")+theme_classic() +
  scale_fill_manual(values = c("#08415c", "#cc2936"), breaks = c("nCRT-R", "nCRT-NR")) +
  geom_hline(yintercept= median(sampleTableCount.Clinical$MATH_score), linetype="dashed", color = "black")+
  annotate(geom="text", x=50, y=50, label=paste("median MATH= ",sprintf("%.2f",median(sampleTableCount.Clinical$MATH_score))), color="black")+
  xlab("SAMPLES") +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gridExtra::grid.arrange(p, bottom="MATH-score distribution from 69 samples")

aux<-(unique(somatic.69[,c("index", "ExonicFunc.refGene")]))  
aux<-as.data.frame(table(aux$ExonicFunc.refGene))
aux$perc <- (aux$Freq/sum(aux$Freq)) *100
aux

# PLOT BOXPLOT MATH_score CORRELATION CLINICAL
#==============================================================================#
# RESPONSE 1
table(sampleTableCount.Clinical$Response1)
aux<-(table(sampleTableCount.Clinical$Response1))
p <- ggplot(sampleTableCount.Clinical  , aes(x=Response1, y=MATH_score, fill=Response1)) + 
  geom_boxplot(outlier.shape = NA)+ 
  scale_fill_manual(values = c("#08415c", "#cc2936"), breaks = c("nCRT-R", "nCRT-NR")) +
  annotate("text", label = paste('N = ',aux[2], sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = ',aux[1], sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  stat_compare_means()+ theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
table(sampleTableCount.Clinical$Response1)
gridExtra::grid.arrange(p, bottom="MATH-score correlation nCRT-R x nCRT-NR (69 samples)")

# sem ROP-107 e ROP-83
#==============================================================================#
table(sampleTableCount.Clinical2$Response1)
p <- ggplot(sampleTableCount.Clinical2  , aes(x=Response1, y=MATH_score, fill=Response1)) + 
  geom_boxplot(outlier.shape = NA)+ 
  scale_fill_manual(values = c("#08415c", "#cc2936"), breaks = c("nCRT-R", "nCRT-NR")) +
  annotate("text", label = paste('N = 14', sep = ''), x = 2, y = 3, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 53', sep = ''), x = 1, y = 3, angle = 0, size = 3.75) +
  stat_compare_means()+ theme_classic()+
  theme(legend.position="none") + geom_jitter(shape=16, position=position_jitter(0.2))

gridExtra::grid.arrange(p, bottom="MATH-score correlation nCRT-R x nCRT-NR (67 samples)")
table(sampleTableCount.Clinical2$Response1)

my_comparisons <- list( c("nCRT-R", "nCRT-NR Not-Metastatic"), 
                        c("nCRT-NR Not-Metastatic", "nCRT-NR Metastatic"), 
                        c("nCRT-R", "nCRT-NR Metastatic") )
# RESPONSE 2
p <- ggplot(sampleTableCount.Clinical, aes(x=Response2, y=MATH_score, fill=Response2)) + 
  geom_boxplot(outlier.shape = NA)+
  stat_compare_means(label.y = 80)+
  scale_fill_manual(values = c("#08415c","#cc2956", "#ef767a"), breaks = c("nCRT-R", "nCRT-NR Metastatic", "nCRT-NR Not-Metastatic")) +
  annotate("text", label = paste('N = 16', sep = ''), x = 3, y = 3, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 37', sep = ''), x = 2, y = 3, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 16', sep = ''), x = 1, y = 3, angle = 0, size = 3.75) +
  stat_compare_means(comparisons = my_comparisons)+
  theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="MATH-score correlation Response2 (69 samples)")
table(sampleTableCount.Clinical$Response2)

# sem ROP-107 e ROP-83
#==============================================================================#
p <- ggplot(sampleTableCount.Clinical2, aes(x=Response2, y=MATH_score, fill=Response2)) + 
  geom_boxplot(outlier.shape = NA)+
  stat_compare_means(label.y = 80)+
  scale_fill_manual(values = c("#08415c","#cc2956", "#ef767a"), breaks = c("nCRT-R", "nCRT-NR Metastatic", "nCRT-NR Not-Metastatic")) +
  annotate("text", label = paste('N = 16', sep = ''), x = 3, y = 3, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 37', sep = ''), x = 2, y = 3, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 14', sep = ''), x = 1, y = 3, angle = 0, size = 3.75) +
  stat_compare_means(comparisons = my_comparisons)+
  theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="MATH-score correlation Response2 (67 samples)")
table(sampleTableCount.Clinical2$Response2)
# RESPONSE 3
p <- ggplot(sampleTableCount.Clinical  , aes(x=Response3, y=MATH_score, fill=Response3)) + 
  geom_boxplot(outlier.shape = NA)+ stat_compare_means()+ 
  scale_fill_manual(values = c("#cc2956", "#ef764a"), breaks = c("Metastatic", "Not-Metastatic")) +
  annotate("text", label = paste('N = 53', sep = ''), x = 2, y = 3, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 16', sep = ''), x = 1, y = 3, angle = 0, size = 3.75) +
  theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="MATH-score correlation Response2 (69 samples)")
table(sampleTableCount.Clinical$Response3)
# sem ROP-107 e ROP-83
#==============================================================================#
p <- ggplot(sampleTableCount.Clinical2  , aes(x=Response3, y=MATH_score, fill=Response3)) + 
  geom_boxplot(outlier.shape = NA)+ stat_compare_means()+ 
  scale_fill_manual(values = c("#cc2956", "#ef764a"), breaks = c("Metastatic", "Not-Metastatic")) +
  annotate("text", label = paste('N = 51', sep = ''), x = 2, y = 3, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 16', sep = ''), x = 1, y = 3, angle = 0, size = 3.75) +
  theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="MATH-score correlation Response2 (67 samples)")
table(sampleTableCount.Clinical2$Response3)


## BARPLOT TMB ####
#==============================================================================#
## TMB colour by Resposta1 
p<- ggplot(data=sampleTableCount.Clinical, aes(x= reorder(ID_Exoma, -TMB), y=TMB, fill=Response1)) +
  geom_bar(stat="identity")+theme_classic() +
  scale_y_sqrt()+
  scale_fill_manual(values = c("#08415c", "#cc2936"), breaks = c("nCRT-R", "nCRT-NR")) +
  geom_hline(yintercept= median(sampleTableCount.Clinical$TMB), linetype="dashed", color = "black")+
  annotate(geom="text", x=50, y=50, label=paste("median TMB= ",sprintf("%.2f",median(sampleTableCount.Clinical$TMB))), color="black")+
  xlab("SAMPLES")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gridExtra::grid.arrange(p, bottom="TMB distribution from 69 samples")

sampleTableCount.Clinical$Perc_TMB <- sampleTableCount.Clinical$TMB/max(sampleTableCount.Clinical$TMB)

head(sampleTableCount.Clinical[,c(8,4,7,40,5,37,38,39)])

# PLOT BOXPLOT TMB CORRELATION CLINICAL
#==============================================================================#
# RESPONSE 1
p <- ggplot(sampleTableCount.Clinical  , aes(x=Response1, y=TMB, fill=Response1)) + 
  geom_boxplot(outlier.shape = NA)+ 
  scale_y_sqrt()+
  scale_fill_manual(values = c("#08415c", "#cc2936"), breaks = c("nCRT-R", "nCRT-NR")) +
  annotate("text", label = paste('N = 16', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 53', sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  stat_compare_means()+ theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="TMB correlation Response1 (69 samples)")
table(sampleTableCount.Clinical$Response1)
# sem ROP-107 e ROP-83
#==============================================================================#
p <- ggplot(sampleTableCount.Clinical2  , aes(x=Response1, y=TMB, fill=Response1)) + 
  geom_boxplot(outlier.shape = NA)+ 
  # scale_y_sqrt()+
  scale_fill_manual(values = c("#08415c", "#cc2936"), breaks = c("nCRT-R", "nCRT-NR")) +
  annotate("text", label = paste('N = 14', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 53', sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  stat_compare_means()+ theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="TMB correlation Response1 (67 samples)")
table(sampleTableCount.Clinical2$Response1)

my_comparisons <- list( c("nCRT-R", "nCRT-NR Not-Metastatic"), 
                        c("nCRT-NR Not-Metastatic", "nCRT-NR Metastatic"), 
                        c("nCRT-R", "nCRT-NR Metastatic") )
# RESPONSE 2
p <- ggplot(sampleTableCount.Clinical, aes(x=Response2, y=TMB, fill=Response2)) + 
  geom_boxplot(outlier.shape = NA)+
  scale_y_sqrt()+
  stat_compare_means(label.y = 10)+
  stat_compare_means(comparisons = my_comparisons)+ 
  scale_fill_manual(values = c("#08415c","#cc2956", "#ef767a"), breaks = c("nCRT-R", "nCRT-NR Metastatic", "nCRT-NR Not-Metastatic")) +
  annotate("text", label = paste('N = 16', sep = ''), x = 3, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 37', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 16', sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="TMB correlation Response2 (69 samples)")
table(sampleTableCount.Clinical$Response2)
# sem ROP-107 e ROP-83
#==============================================================================#
p <- ggplot(sampleTableCount.Clinical2, aes(x=Response2, y=TMB, fill=Response2)) + 
  geom_boxplot(outlier.shape = NA)+
  # scale_y_sqrt()+
  stat_compare_means(label.y = 10)+
  stat_compare_means(comparisons = my_comparisons)+ 
  scale_fill_manual(values = c("#08415c","#cc2956", "#ef767a"), breaks = c("nCRT-R", "nCRT-NR Metastatic", "nCRT-NR Not-Metastatic")) +
  annotate("text", label = paste('N = 14', sep = ''), x = 3, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 37', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 16', sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="TMB correlation Response2 (67 samples)")
table(sampleTableCount.Clinical2$Response2)

# RESPONSE 3
p <- ggplot(sampleTableCount.Clinical, aes(x=Response3, y=TMB, fill=Response3)) + 
  geom_boxplot(outlier.shape = NA)+
  scale_y_sqrt()+
  stat_compare_means(label.y = 10)+
  scale_fill_manual(values = c("#cc2956", "#ef764a"), breaks = c("Metastatic", "Not-Metastatic")) +
  annotate("text", label = paste('N = 53', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 16', sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="TMB correlation Response3 (69 samples)")
table(sampleTableCount.Clinical$Response3)
# sem ROP-107 e ROP-83
#==============================================================================#
p <- ggplot(sampleTableCount.Clinical2, aes(x=Response3, y=TMB, fill=Response3)) + 
  geom_boxplot(outlier.shape = NA)+
  # scale_y_sqrt()+
  stat_compare_means(label.y = 10)+
  scale_fill_manual(values = c("#cc2956", "#ef764a"), breaks = c("Metastatic", "Not-Metastatic")) +
  annotate("text", label = paste('N = 51', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 16', sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="TMB correlation Response3 (67 samples)")
table(sampleTableCount.Clinical2$Response3)


#==============================================================================#
#                 CORRELACAO TMB vs MATH-SCORE           ####
#==============================================================================#

shapiro.test(sampleTableCount.Clinical$TMB)
qqnorm(sampleTableCount.Clinical$TMB)
qqline(sampleTableCount.Clinical$TMB)

shapiro.test(sampleTableCount.Clinical$MATH_score)
qqnorm(sampleTableCount.Clinical$MATH_score)
qqline(sampleTableCount.Clinical$MATH_score)


sampleTableCount.Clinical2 <- subset(sampleTableCount.Clinical, !(ID_Exoma %in% c("ROP-107", "ROP-83")), )

shapiro.test(sampleTableCount.Clinical2$TMB)
qqnorm(sampleTableCount.Clinical2$TMB)
qqline(sampleTableCount.Clinical2$TMB)

shapiro.test(sampleTableCount.Clinical2$MATH_score)
qqnorm(sampleTableCount.Clinical2$MATH_score)
qqline(sampleTableCount.Clinical2$MATH_score)

p<-ggscatter(sampleTableCount.Clinical, x = "MATH_score", y = "TMB",
             add = "reg.line",  # Add regressin line
             # add = "loess",  # Local Polynomial Regression Fitting
             #          cor.coef=T,
             cor.method = "spearman",
             add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
             conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",  label.x = 30,label.y = 100)
gridExtra::grid.arrange(p, bottom="Spearman Test (69 samples)")

sp <- ggscatter(sampleTableCount.Clinical2, x = "MATH_score", y = "TMB",
                add = "reg.line",  # Add regressin line
                # add = "loess",  # Local Polynomial Regression Fitting
                # cor.coef=T,
                cor.method = "spearman",
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 30, label.y = 17)
gridExtra::grid.arrange(sp, bottom="Spearman Test (67 samples)")

#sem Metastatic
sp <- ggscatter(subset(sampleTableCount.Clinical2,Response3 == "Not-Metastatic"), x = "MATH_score", y = "TMB",
                add = "reg.line",  # Add regressin line
                cor.method = "spearman",
                color = "white",
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 40)+ 
  geom_point(aes(color = Response2))+
  theme(legend.position="top")
gridExtra::grid.arrange(sp, bottom="Spearman Test (67 samples) excluded Metastatic")


# RESPONSE 1

sp <- ggscatter(sampleTableCount.Clinical2, x = "MATH_score", y = "TMB",
                add = "reg.line",  # Add regressin line
                cor.method = "spearman",
                color = "white",
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman",  aes(color = Response1), label.x = 40)+
  geom_point(aes(color = Response1))+
  theme(legend.position="top")
gridExtra::grid.arrange(sp, bottom="Spearman Test (67 samples)")

# RESPONSE 2
sp <- ggscatter(sampleTableCount.Clinical2, x = "MATH_score", y = "TMB",
                add = "reg.line",  # Add regressin line
                cor.method = "spearman",
                color = "white",
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman",  aes(color = Response2), label.x = 40)+
  geom_point(aes(color = Response2))+
  theme(legend.position="top")
gridExtra::grid.arrange(sp, bottom="Spearman Test (67 samples)")

# RESPONSE 2 NOVO
sp <- ggscatter(sampleTableCount.Clinical2, x = "MATH_score", y = "TMB",
                add = "reg.line",  # Add regressin line
                cor.method = "spearman",
                color = "white",
                add.params = list(color = "Response2", fill = "Response2"), # Customize reg. line
                conf.int = F # Add confidence interval
) + stat_cor(method = "spearman",  aes(color = Response2,fill = "Response2"), label.x = 40)+
  geom_point(aes(color = Response2))+
  theme(legend.position="top")
gridExtra::grid.arrange(sp, bottom="Spearman Test (67 samples)")

# RESPONSE 2  excluded metastatic
sp <- ggscatter(subset(sampleTableCount.Clinical2,Response2 != "nCRT-NR Metastatic"), x = "MATH_score", y = "TMB",
                add = "reg.line",  # Add regressin line
                cor.method = "spearman",
                color = "white",
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman",  aes(color = Response2), label.x = 40)+
  geom_point(aes(color = Response2))+
  theme(legend.position="top")
gridExtra::grid.arrange(sp, bottom="Spearman Test (67 samples) excluded Metastatic")
# RESPONSE 3
sp <- ggscatter(sampleTableCount.Clinical2, x = "MATH_score", y = "TMB",
                add = "reg.line",  # Add regressin line
                cor.method = "spearman",
                color = "white",
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman",  aes(color = Response3), label.x = 40)+
  geom_point(aes(color = Response3))+
  theme(legend.position="top")
gridExtra::grid.arrange(sp, bottom="Spearman Test (67 samples)")


#==============================================================================#
# TOP 50 GENES MAIS MUTADOS ####
#==============================================================================#
dgenes<-fread("input_Somatic-Filter/driver-genes.txt", header = F, col.names = "gene")
dMMR<-fread("input_Somatic-Filter/dMMR_genes.txt", header = F, col.names = "gene")

mut.genes=unique(select(somatic.69, c("index","Gene.refGene")))
aux2 <- as.data.frame(table(mut.genes$Gene.refGene))
names(aux2) <- c("Gene", "Freq_Mutation")
aux2 <- (aux2[order(aux2$Freq, decreasing = T),])
aux2$class <-'none'
aux2$class[aux2$Gene %in% dgenes$gene] <-'driver'
aux2$class[aux2$Gene %in% dMMR$gene] <-'dMMR'
#fwrite(aux2, "output_Somatic-Filter/countMut_by_Genes.tsv", quote = F, sep="\t")
aux2=head(aux2,50)
top50<- aux2$Gene
p<-ggplot(data=aux2, aes(x= reorder(Gene, -Freq_Mutation), y=Freq_Mutation, fill=class)) +
  geom_bar(stat="identity" )+
  # scale_y_sqrt()+
  theme_classic() +xlab("Gene")
p +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#==============================================================================#
# TOP 50 GENES MAIS MUTADOS 67 samples ####
#==============================================================================#
dgenes<-fread("input_Somatic-Filter/driver-genes.txt", header = F, col.names = "gene")
dMMR<-fread("input_Somatic-Filter/dMMR_genes.txt", header = F, col.names = "gene")
somatic.67 <- subset(somatic.69, !(ID_Exoma %in% c("ROP-83", "ROP-107")), )

mut.genes=unique(select(somatic.67, c("index","Gene.refGene")))
aux2 <- as.data.frame(table(mut.genes$Gene.refGene))
names(aux2) <- c("Gene", "Freq_Mutation")
aux2 <- (aux2[order(aux2$Freq, decreasing = T),])
aux2$class <-'none'
aux2$class[aux2$Gene %in% dgenes$gene] <-'driver'
aux2$class[aux2$Gene %in% dMMR$gene] <-'dMMR'
#fwrite(aux2, "output_Somatic-Filter/countMut_by_Genes.tsv", quote = F, sep="\t")
aux2=head(aux2,50)
top50<- aux2$Gene
p<-ggplot(data=aux2, aes(x= reorder(Gene, -Freq_Mutation), y=Freq_Mutation, fill=class)) +
  geom_bar(stat="identity" )+
  # scale_y_sqrt()+
  theme_classic() +xlab("Gene")
p +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#==============================================================================#
# DNDS_CV ####

library(dndscv)

#==============================================================================#
# TESTE PARA CORRIGIR A ANALISE DE dNdS_cv PELA REGIÃO DO EXOMA
#==============================================================================#

# path_genome_fasta = system.file("extdata", "chr3_segment.fa", package = "dndscv", mustWork = TRUE)
# 
# target_file<- fread("input_Somatic-Filter/xgen-exome-research-panel-v2-targets-hg38.bed ")
# target_file$Gene_v <-target_file$V4
# target_file <- separate(data = target_file, col = V4, into = c("V4"), sep = "_")
# target_file <- target_file[,c(4,4,7,1,2,3,2,3,3-2,6)]
# colnames(target_file) <- c("gene.id", "gene.name", "cds.id", "chr", "chr.coding.start", "chr.coding.end", "cds.start", "cds.end", "length", "strand")
# target_file[target_file$strand == "+",]$strand <- "1"
# target_file[target_file$strand == "-",]$strand <- "-1"
# fwrite(target_file, "input_Somatic-Filter/xgen-exome-research-panel-v2-targets-hg38.cds", quote = F, col.names = T, row.names = F)
# fwrite(as.data.frame(unique(target_file$gene.name)), "input_Somatic-Filter/xgen-exome-research-panel-v2-targets-hg38.list", quote = F, col.names = F, row.names = F)
# 

# cd input_Somatic-Filter/
# cut -f 4 xgen-exome-research-panel-v2-targets-hg38.bed | cut -d '_' -f1| sort -u > xgen-exome-research-panel-v2-targets-hg38.list
# path_cds_table <- c("input_Somatic-Filter/xgen-exome-research-panel-v2-targets-hg38.mart_export.txt")
# path_cds_table <- c("input_Somatic-Filter/mart_export.txt")
# path_cds_table <- c("input_Somatic-Filter/chr19-hg38.mart_export.txt")
# 
# path_genome_fasta = c("input_Somatic-Filter/Homo_sapiens.GRCh38.dna.chromosome.19.fa")
# path_genome_fasta = c("input_Somatic-Filter/Homo_sapiens.GRCh38.p14_GCA_000001405.29.fa")
# buildref(cdsfile=path_cds_table, genomefile=path_genome_fasta, outfile = "example_output_refcds.rda", excludechrs="MT", useids = T,)

#==============================================================================#
# 69 SAMPLES
#==============================================================================#
all.mut.flt<-fread(paste0(output,"/mutation_codRegion_Freq005.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)
all.mut.flt.69<- subset(all.mut.flt, SAMPLE %in% unique(Clinical$Sample))
# 
# aux<-(mutations2[which(mutations2$Ref %in% c("A", "T","C", "G")),])
# aux<-(aux[which(aux$Alt %in% c("A", "T","C", "G")),])
#  unique(aux[,c('Ref', 'Alt')])

# targert_file<- fread("input_Somatic-Filter/xgen-exome-research-panel-v2-targets-hg38.bed ")
# targert_file <- separate(data = targert_file, col = V4, into = c("Gene","id"), sep = "_")
mutations <- all.mut.flt.69[,c("SAMPLE", "Chr", "Start", "Ref", "Alt")]
dndsout = dndscv(mutations, refdb="../resultados_Mutect2/input_Somatic-Filter/RefCDS_human_GRCh38.p12_dNdScv.0.1.0.rda", cv=NULL,
                 # gene_list = unique(targert_file$Gene),
                 max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf )
mutations2 <-mutations
mutations2$Chr <- gsub("chr", "", mutations2$Chr)
dndsout2 = dndscv(mutations2, refdb="../resultados_Mutect2/input_Somatic-Filter/RefCDS_human_GRCh38_GencodeV18_recommended.rda", cv=NULL,
                  # gene_list = unique(targert_file$Gene),
                  max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf )

mutations3 <- all.mut.flt.69[,c("SAMPLE", "CHR", "POS", "REF", "ALT")]
dndsout3 = dndscv(mutations3, refdb="../resultados_Mutect2/input_Somatic-Filter/RefCDS_human_GRCh38.p12_dNdScv.0.1.0.rda", cv=NULL,
                  # gene_list = unique(targert_file$Gene),
                  max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf )
mutations4 <-mutations3
mutations4$CHR <- gsub("chr", "", mutations4$CHR)
dndsout4 = dndscv(mutations2, refdb="../resultados_Mutect2/input_Somatic-Filter/RefCDS_human_GRCh38_GencodeV18_recommended.rda", cv=NULL,
                  # gene_list = unique(targert_file$Gene),
                  max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf )


for (i in list(dndsout, dndsout2, dndsout3, dndsout4)) {
  print(i$globaldnds)
  print(i$nbreg$theta)
  
  signif_genes_localmodel = as.vector(i$sel_loc$gene_name[i$sel_loc$qall_loc<0.01])
  print(signif_genes_localmodel)
  # View(dndsout2$sel_cv) # This is shown as an example but these results based on a few genes should not be trusted
  sel_cv = i$sel_cv
  signif_genes = sel_cv[sel_cv$pglobal_cv<=0.05, c("gene_name", "pglobal_cv" ,"qglobal_cv")]
  #print(subset(signif_genes,qglobal_cv <= 0.05))
  print(nrow(subset(signif_genes,pglobal_cv <= 0.05)))
  print(nrow(subset(signif_genes,qglobal_cv <= 0.1)))
  print(nrow(subset(signif_genes,qglobal_cv <= 0.05)))
  
}

print(dndsout2$globaldnds)
print(dndsout2$nbreg$theta)

signif_genes_localmodel = as.vector(dndsout2$sel_loc$gene_name[dndsout2$sel_loc$qall_loc<0.01])
print(signif_genes_localmodel)
# View(dndsout2$sel_cv) # This is shown as an example but these results based on a few genes should not be trusted
sel_cv = dndsout2$sel_cv
signif_genes = sel_cv[sel_cv$pglobal_cv<=0.05, c("gene_name", "pglobal_cv" ,"qglobal_cv")]
#print(subset(signif_genes,qglobal_cv <= 0.05))
nrow(subset(signif_genes,pglobal_cv <= 0.05))
nrow(subset(signif_genes,qglobal_cv <= 0.1))
nrow(subset(signif_genes,qglobal_cv <= 0.05))

fwrite(sel_cv[sel_cv$pglobal_cv<0.05, ], paste0(output,"/somatic.69.mut_signif_dndscv.tsv"), quote = F, sep="\t")
fwrite(signif_genes, paste0(output,"/somatic.69.mut_signif_genes.tsv"), quote = F, sep="\t")

#==============================================================================#
#                ONCOPLOT  69 samples ####
#==============================================================================#
signif_genes <- fread(file=paste0(output,"/somatic.69.mut_signif_genes.tsv"))

somatic.maf<-fread(paste0(output,"/mutation_codRegion_Freq005.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)
somatic.maf.69<- subset(somatic.maf, SAMPLE %in% unique(Clinical$Sample))
somatic.maf.69 <- merge(somatic.maf.69, sampleTableCount.Clinical[,c("SAMPLE","ID_Exoma")], by= "SAMPLE", all.x = T)
# somatic.maf<- subset(somatic.maf.69, select=c("CHR", "POS", "REF", "ALT", "Gene.refGene","ID_Exoma","ExonicFunc.refGene"))
somatic.maf<- subset(somatic.maf.69, select=c("Chr", "Start",  "End","Ref", "Alt", "Gene.refGene","ID_Exoma","ExonicFunc.refGene"))

somatic.maf<- somatic.maf %>% 
  mutate(Gene.refGene = strsplit(as.character(Gene.refGene), ";")) %>% 
  unnest(Gene.refGene)

# DESCOMENTE ABAIXO PARA PLOTAR SEM MUTAÇOES SINONIMAS
# somatic.maf<- subset(somatic.69, select=c("Chr", "Start",  "End","Ref", "Alt", "Gene.refGene","ID_Exoma","ExonicFunc.refGene"))
# somatic.maf<- somatic.maf %>%
#   mutate(Gene.refGene = strsplit(as.character(Gene.refGene), ";")) %>%
#   unnest(Gene.refGene)

somatic.maf$Reference_Allele <- somatic.maf$Ref
colnames(somatic.maf)[names(somatic.maf)=="Ref"]<-"Tumor_Seq_Allele1"
colnames(somatic.maf)[names(somatic.maf)=="Alt"]<-"Tumor_Seq_Allele2"
colnames(somatic.maf)[names(somatic.maf)=="Gene.refGene"]<-"Hugo_Symbol"
colnames(somatic.maf)[names(somatic.maf)=="Chr"]<-"Chromosome"
colnames(somatic.maf)[names(somatic.maf)=="Start"]<-"Start_Position"
colnames(somatic.maf)[names(somatic.maf)=="End"]<-"End_Position"
colnames(somatic.maf)[names(somatic.maf)=="ExonicFunc.refGene"]<-"Variant_Classification"
somatic.maf$Variant_Type<-"ExonicFunc.refGene"
colnames(somatic.maf)[names(somatic.maf)=="ID_Exoma"]<-"Tumor_Sample_Barcode"

sampleTableCount.Clinical.maf <- sampleTableCount.Clinical
colnames(sampleTableCount.Clinical.maf)[names(sampleTableCount.Clinical.maf)=="ID_Exoma"]<-"Tumor_Sample_Barcode"

aux<-unique(subset(somatic.69, class == "dMMR",select = c("SAMPLE","class") ))
sampleTableCount.Clinical.maf$dMMR = "pMMR"
sampleTableCount.Clinical.maf$dMMR[sampleTableCount.Clinical.maf$SAMPLE %in% aux$SAMPLE] <-'dMMR'

maf = read.maf(maf = somatic.maf,
               clinicalData = sampleTableCount.Clinical.maf,
               vc_nonSyn = unique(somatic.maf$Variant_Classification),
               verbose = F)
unique(somatic.maf$Variant_Classification)

oncoplot(maf = maf, 
         genes = subset(signif_genes, qglobal_cv <= 0.1, select=gene_name)$gene_name,
         sortByAnnotation = T,
         clinicalFeatures = c("Response1","Response2","Response3"), 
         showTumorSampleBarcodes = T, SampleNamefontSize = 0.8,
         removeNonMutated = F,  drawColBar = T,
         fontSize = 0.6, annoBorderCol = "white",
         annotationFontSize = 1, 
         legend_height =2, 
         legendFontSize = 1)


#==============================================================================#
# DNDS_CV 2 - sem ROP-107 e ROP-83####
#==============================================================================#
all.mut.flt<-fread(paste0(output,"/mutation_codRegion_Freq005.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)
all.mut.flt.69<- subset(all.mut.flt, SAMPLE %in% unique(Clinical$Sample))
all.mut.flt.67<-subset(all.mut.flt.69, 
                       !(SAMPLE %in% c("ROP-83-ExC85-xgenV2_S54", "ROP-107-ExC85-xgenV2_S71")),)
mutations2 <- all.mut.flt.67[,c("SAMPLE", "Chr", "Start", "Ref", "Alt")]
mutations2$Chr <- gsub("chr", "", mutations2$Chr)
dndsout2 = dndscv(mutations2, refdb="../resultados_Mutect2/input_Somatic-Filter/RefCDS_human_GRCh38_GencodeV18_recommended.rda", cv=NULL,
                  # gene_list = unique(targert_file$Gene),
                  max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf )

print(dndsout2$globaldnds)
print(dndsout2$nbreg$theta)
signif_genes_localmodel = as.vector(dndsout$sel_loc$gene_name[dndsout$sel_loc$qall_loc<0.1])
print(signif_genes_localmodel)
# View(dndsout$sel_cv) # This is shown as an example but these results based on a few genes should not be trusted
sel_cv = dndsout2$sel_cv
signif_genes = sel_cv[sel_cv$pglobal_cv<0.05, c("gene_name", "pglobal_cv" ,"qglobal_cv")]
print(subset(signif_genes,qglobal_cv <= 0.1))
nrow(subset(signif_genes,pglobal_cv <= 0.05))
nrow(subset(signif_genes,qglobal_cv <= 0.1))
nrow(subset(signif_genes,qglobal_cv <= 0.05))
fwrite(sel_cv[sel_cv$pglobal_cv<0.05, ], paste0(output,"/somatic.67.mut_signif_dndscv.tsv"), quote = F, sep="\t")
fwrite(signif_genes, paste0(output,"/somatic.67.mut_signif_genes.tsv"), quote = F, sep="\t")

#==============================================================================#
#                ONCOPLOT  67 samples ####
#==============================================================================#
signif_genes <- fread(file=paste0(output,"/somatic.67.mut_signif_genes.tsv"))

somatic.maf<- subset(somatic.67, select=c("CHR", "POS", "REF", "ALT", "Gene.refGene","ID_Exoma","ExonicFunc.refGene"))

somatic.maf<-fread(paste0(output,"/mutation_codRegion_Freq005.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)
somatic.maf.69<- subset(somatic.maf, SAMPLE %in% unique(Clinical$Sample))
somatic.maf.67<- subset(somatic.maf.69, 
                        !(SAMPLE %in% c("ROP-83-ExC85-xgenV2_S54", "ROP-107-ExC85-xgenV2_S71")),)
somatic.maf.67 <- merge(somatic.maf.67, sampleTableCount.Clinical2[,c("SAMPLE","ID_Exoma")], by= "SAMPLE", all.x = T)
somatic.maf<- subset(somatic.maf.67, select=c("Chr", "Start",  "End","Ref", "Alt", "Gene.refGene","ID_Exoma","ExonicFunc.refGene"))

somatic.maf<- somatic.maf %>% 
  mutate(Gene.refGene = strsplit(as.character(Gene.refGene), ";")) %>% 
  unnest(Gene.refGene)

# DESCOMENTE ABAIXO PARA PLOTAR SEM MUTAÇOES SINONIMAS
# somatic.maf<- subset(somatic.67, select=c("Chr", "Start",  "End","Ref", "Alt", "Gene.refGene","ID_Exoma","ExonicFunc.refGene"))
# somatic.maf<- somatic.maf %>%
#   mutate(Gene.refGene = strsplit(as.character(Gene.refGene), ";")) %>%
#   unnest(Gene.refGene)

somatic.maf$Reference_Allele <- somatic.maf$Ref
colnames(somatic.maf)[names(somatic.maf)=="Ref"]<-"Tumor_Seq_Allele1"
colnames(somatic.maf)[names(somatic.maf)=="Alt"]<-"Tumor_Seq_Allele2"
colnames(somatic.maf)[names(somatic.maf)=="Gene.refGene"]<-"Hugo_Symbol"
colnames(somatic.maf)[names(somatic.maf)=="Chr"]<-"Chromosome"
colnames(somatic.maf)[names(somatic.maf)=="Start"]<-"Start_Position"
colnames(somatic.maf)[names(somatic.maf)=="End"]<-"End_Position"
colnames(somatic.maf)[names(somatic.maf)=="ExonicFunc.refGene"]<-"Variant_Classification"
somatic.maf$Variant_Type<-"ExonicFunc.refGene"
colnames(somatic.maf)[names(somatic.maf)=="ID_Exoma"]<-"Tumor_Sample_Barcode"

sampleTableCount.Clinical.maf <- sampleTableCount.Clinical2
colnames(sampleTableCount.Clinical.maf)[names(sampleTableCount.Clinical.maf)=="ID_Exoma"]<-"Tumor_Sample_Barcode"

aux<-unique(subset(somatic.67, class == "dMMR",select = c("SAMPLE","class") ))
sampleTableCount.Clinical.maf$dMMR = "pMMR"
sampleTableCount.Clinical.maf$dMMR[sampleTableCount.Clinical.maf$SAMPLE %in% aux$SAMPLE] <-'dMMR'

maf = read.maf(maf = somatic.maf,
               clinicalData = sampleTableCount.Clinical.maf,
               vc_nonSyn = unique(somatic.maf$Variant_Classification),
               verbose = F)
unique(somatic.maf$Variant_Classification)

oncoplot(maf = maf, 
         genes = subset(signif_genes, qglobal_cv <= 0.1, select=gene_name)$gene_name,
         sortByAnnotation = T,
         clinicalFeatures = c("Response1","Response2","Response3"), 
         showTumorSampleBarcodes = T, SampleNamefontSize = 0.8,
         removeNonMutated = F,  drawColBar = T,
         fontSize = 0.6, annoBorderCol = "white",
         annotationFontSize = 1, 
         legend_height =2, 
         legendFontSize = 1)


#==============================================================================#
#  pheatmap - GENES with dMMR  
#==============================================================================#
mut.genes=unique(select(somatic.69, c("Gene.refGene","ID_Exoma", "class")))
aux<- as.data.frame(unique(subset(mut.genes, select= ID_Exoma)))

mut.genes.heatmap= as.data.frame(table(unique(select(somatic.69, c("Gene.refGene","ID_Exoma")))))
mut.genes.heatmap<- mut.genes.heatmap[mut.genes.heatmap$Gene.refGene %in% dMMR$gene,]
mut.genes.heatmap<- mut.genes.heatmap[mut.genes.heatmap$ID_Exoma %in% aux$ID_Exoma,]
df<- as.matrix(reshape2::dcast(mut.genes.heatmap, Gene.refGene ~ ID_Exoma,  value.var = "Freq"))
rownames(df) <- df[,1]
df<-type.convert(df[,-1], as.is = T)

## Juntando tabelas TMB_MATH_#MUT +  dados_clinicos 
my_sample_col <- merge(sampleTableCount, Clinical, by.x = "SAMPLE", by.y = "Sample")
my_sample_col <- data.frame(subset(my_sample_col, select = c("ID_Exoma","Response1","Response2","Response3")))

row.names(my_sample_col) <- my_sample_col[,1]
my_sample_col<-subset(my_sample_col, select=-1)

my_sample_col<-my_sample_col[order(my_sample_col$Response1 ),]
df_ordered <- df[,rownames(my_sample_col)]

pheatmap((df_ordered), color=colorRampPalette(c("white", "darkred"))(100),
         annotation_col = my_sample_col,
         border_color = "grey60",
         cluster_rows  = F,cluster_cols  = F, fontsize = 8,
         legend = F )


#==============================================================================#
#  pheatmap - GENES with dMMR  67samples
#==============================================================================#
# mut.genes=unique(select(somatic.67, c("Gene.refGene","ID_Exoma", "class")))
# aux<- as.data.frame(unique(subset(mut.genes, select= ID_Exoma)))
# 
# mut.genes.heatmap= as.data.frame(table(unique(select(somatic.67, c("Gene.refGene","ID_Exoma")))))
# mut.genes.heatmap<- mut.genes.heatmap[mut.genes.heatmap$Gene.refGene %in% dMMR$gene,]
# mut.genes.heatmap<- mut.genes.heatmap[mut.genes.heatmap$ID_Exoma %in% aux$ID_Exoma,]
# df<- as.matrix(reshape2::dcast(mut.genes.heatmap, Gene.refGene ~ ID_Exoma,  value.var = "Freq"))
# rownames(df) <- df[,1]
# df<-type.convert(df[,-1], as.is = T)
# 
# my_sample_col <- sampleTableCount.Clinical2
# 
# my_sample_col <- data.frame(subset(my_sample_col, select = c("ID_Exoma","Response1","Response2","Response3")))
# row.names(my_sample_col) <- my_sample_col[,1]
# my_sample_col<-subset(my_sample_col, select=-1)
# 
# my_sample_col<-my_sample_col[order(my_sample_col$Response1 ),]
# df_ordered <- df[,rownames(my_sample_col)]
# 
# pheatmap((df_ordered), color=colorRampPalette(c("white", "darkred"))(100),
#          annotation_col = my_sample_col,
#          border_color = "grey60",
#          cluster_rows  = F,cluster_cols  = F, fontsize = 8,
#          legend = F )
# 







#==============================================================================#
# PIPELINE EXECUTADO CORRETAMENTE ATE AQUI,
# POIS EU NAO RODEI A ANALISE DE ASSINATURA MUTACIONAL (MUTALISK)
#==============================================================================#







#==============================================================================#
#        MUTALISK - ASSINATURAS MUTACIONAIS      ####
#==============================================================================#
#RODAR MUTALISK E BAIXAR PARA DENTRO DO DIRETORIO

# comando SHELL para cria a tabela
# grep -e '**THIS_SIGs\|**THIS_SIG_CONTRIBUTIONS\|**Cosine_Similarity' Mutalisk/*.txt | sed -E 's/_user.select_[0-9]_report.txt\:\*\*/\t/' | sed 's/Mutalisk\///' | sed 's/Cosine_Similarity\t/Cosine_Similarity/'| sed 's/\t/ /' |cut -d " " -f 1,2 > table_signatures.1.tsv
# grep -e '**THIS_SIGs\|**THIS_SIG_CONTRIBUTIONS\|**Cosine_Similarity' Mutalisk/*.txt | sed -E 's/_user.select_[0-9]_report.txt\:\*\*/\t/' | sed 's/Mutalisk\///' | sed 's/Cosine_Similarity\t/Cosine_Similarity/'| sed 's/\t/ /' |cut -d " " -f 3- | sed 's/ /;/g' > table_signatures.2.tsv

# grep -e '**THIS_SIGs\|**THIS_SIG_CONTRIBUTIONS\|**Cosine_Similarity' Mutalisk_synonymous/*.txt | sed -E 's/_user.selection_[0-9]_report.txt\:\*\*/\t/' | sed 's/Mutalisk_synonymous\///' | sed 's/Cosine_Similarity\t/Cosine_Similarity/'| sed 's/\t/ /' |cut -d " " -f 1,2 > table_signatures.1.synonymous.tsv
# grep -e '**THIS_SIGs\|**THIS_SIG_CONTRIBUTIONS\|**Cosine_Similarity' Mutalisk_synonymous/*.txt | sed -E 's/_user.selection_[0-9]_report.txt\:\*\*/\t/' | sed 's/Mutalisk_synonymous\///' | sed 's/Cosine_Similarity\t/Cosine_Similarity/'| sed 's/\t/ /' |cut -d " " -f 3- | sed 's/ /;/g' > table_signatures.2.synonymous.tsv


tab1<-fread(paste0(output,"table_signatures.1.tsv"), header = F)
tab2<-fread(paste0(output,"table_signatures.2.tsv"), sep="\t", header = F)

assinaturas<- cbind(tab1, tab2)
colnames(assinaturas) <- c("SAMPLE", "mutalisk", "value")
assinaturas<-dcast(assinaturas, SAMPLE ~ mutalisk)
assinaturas$SAMPLE <- gsub(".recode.final", "", assinaturas$SAMPLE)

ass1<- assinaturas[,c(1,3,4)]
final_Sig<-data.frame()
contador<-1
for (x in 1:nrow(ass1)) {
  sig_id<-as.array(str_split(ass1[x,3], pattern  =  ";"))
  sig_perc<-as.array(str_split(ass1[x,2], pattern  =  ";"))
  for (y in 1:length(sig_id[[1]])) {
    final_Sig[contador,1] <-ass1[x,1]
    final_Sig[contador,2] <-sig_id[[1]][y]
    final_Sig[contador,3] <-sig_perc[[1]][y]
    contador=contador+1
  }
}

colnames(final_Sig)<-c("SAMPLE", "sig_id", "perc")
final_Sig$SAMPLE <- gsub(".recode.final", "", final_Sig$SAMPLE)

final_Sig<-merge(final_Sig, assinaturas[,c(1,2)], by = "SAMPLE" )

final_Sig$dMMR <-"pMMR"
mutados<- unique(subset(somatic.69, class == "dMMR", select = SAMPLE))
final_Sig[final_Sig$SAMPLE %in% mutados$SAMPLE,"dMMR"] <-"dMMR"
# final_Sig<-merge(final_Sig, subset(somatic_Samples.resp, select = c("SAMPLE","Response1","Response3")), by="sample", all.y = T)
final_Sig<-merge(final_Sig, subset(sampleTableCount.Clinical, select = c("SAMPLE","Response1","Response3")), by="SAMPLE", all.y = T)
final_Sig$perc <-as.numeric(final_Sig$perc)
final_Sig$cosmic_sig <- "sig_others"
final_Sig[final_Sig$sig_id %in% c("6","15","20","26"),"cosmic_sig"] <-"MMR_sig"

fwrite(final_Sig, "table_signatures.tsv", sep = "\t", quote = F, col.names = T, row.names = F)

mutsig <-fread("table_signatures.tsv", header = T)
mutsig <- merge(Clinical[,c(1,2)], mutsig, by.x= "Sample",by.y = "SAMPLE")
#mutsig<-subset(mutsig, dMMR == "dMMR")
mutsig$sig_id<-as.character(mutsig$sig_id)

length(unique(mutsig$ID_Exoma))
length(unique(mutsig$Sample))
mutsig <- subset(mutsig, ID_Exoma %in% rownames(my_sample_col))
mutsig$ID_Exoma <- factor(mutsig$ID_Exoma, levels=rownames(my_sample_col))
sort(as.numeric(unique(mutsig$sig_id)))

p<-ggplot(data=mutsig, aes(x=ID_Exoma, y=perc, fill=sig_id)) +
  geom_bar(position="stack", stat="identity") + theme_classic() +xlab("")+
  scale_fill_manual(name="Mutational Signatures Cosmic_v2",
                    #breaks = my_sig_id,
                    breaks = c("10","6","15","20", "3","1","2","4","5","7","8","9","12","13","16","17","18","21","22","24","28","29","30"),
                    values=c("#f3722c","forestgreen","darkgreen","yellowgreen", "#edE2f9","#edf2f4","#fefae0","#edede9","#d6ccc2","#e3d5ca","#fec5bb","#fae1dd","#c9ada7","#ffc999","#ffcdb2","#ffd7ba","#fec89a","#d5bdaf","#cad2c5","#edede9","#d6ccc2","#e3d5ca","#e3d5ca" ))

p +theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=0.2), legend.position = "bottom")



#==============================================================================#
#                 dMMR vs MMR_sig vs MATH_score vs TMB                      ####
#==============================================================================#

mutsig <-fread("table_signatures.tsv", header = T)
mutsig <- merge(Clinical[,c(1,2)], mutsig, by.x= "Sample",by.y = "SAMPLE")

aux<-unique(subset(somatic.69, class == "dMMR",select = c("SAMPLE","class") ))
sampleTableCount.Clinical$dMMR = "pMMR"
sampleTableCount.Clinical$dMMR[sampleTableCount.Clinical$SAMPLE %in% aux$SAMPLE] <-'dMMR'

aux <- as.data.frame(unique(subset(mutsig, cosmic_sig == "MMR_sig", select = "Sample")))
sampleTableCount.Clinical$MMR_Sig = "others_sig"
sampleTableCount.Clinical[sampleTableCount.Clinical$SAMPLE %in% aux$Sample, "MMR_Sig"] <- "MMR_sig"

aux <- sampleTableCount.Clinical

print(table(aux$dMMR))
print(table(aux$MMR_Sig))

print(as.matrix(table(aux[,c("dMMR", "MMR_Sig")]), nr=2))
dt3<- matrix(table(aux[,c("dMMR", "MMR_Sig")]), nr=2)
print(fisher.test(dt3))
print(chisq.test(dt3))

score_MMR_sig <- subset(mutsig, cosmic_sig=="MMR_sig") %>% group_by(Sample, Response1, Response3) %>% 
  dplyr::summarise(score_MMR_sig = sum(perc),
                   .groups = 'drop')

# aux$score_MMR_sig <- 0
aux<- merge(aux, score_MMR_sig[,c(1,4)], by.x= "SAMPLE", by.y= "Sample", all.x = T)
aux$score_MMR_sig[is.na(aux$score_MMR_sig)] <- 0

p<-ggplot(subset(aux, score_MMR_sig>0), aes(x=dMMR, y=score_MMR_sig)) + 
  geom_boxplot(aes(colour = dMMR),outlier.shape = NA) + stat_compare_means() +
  annotate("text", label = paste('N = 25', sep = ''), x = 2, y = -0.05, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 6', sep = ''), x = 1, y = -0.05, angle = 0, size = 3.75) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_classic() +theme(legend.position="none")
gridExtra::grid.arrange(p, bottom="")
print(table(aux[aux$score_MMR_sig >0, ]$dMMR))
print(table(aux$dMMR))

#==============================================================================#
#Response1 vs MMR_sig
p<-ggplot(aux[aux$score_MMR_sig >0, ], 
          aes(x=Response1, y=score_MMR_sig)) + 
  geom_boxplot(outlier.shape = NA) + stat_compare_means() +
  annotate("text", label = paste('N = 10', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 21', sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_classic() +theme(legend.position="left")
gridExtra::grid.arrange(p, bottom="")
print(table(aux[aux$score_MMR_sig >0, ]$Response1))
table(aux$Response1)
#Response2 vs MMR_sig
p<-ggplot(aux[aux$score_MMR_sig >0, ], 
          aes(x=Response2, y=score_MMR_sig)) + 
  geom_boxplot(outlier.shape = NA) + stat_compare_means() +
  annotate("text", label = paste('N = 10', sep = ''), x = 3, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 17', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 4', sep = ''),  x = 1, y = 0, angle = 0, size = 3.75) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  stat_compare_means(comparisons = my_comparisons)+
  theme_classic() +theme(legend.position="left")
gridExtra::grid.arrange(p, bottom="")
print(table(aux[aux$score_MMR_sig >0, ]$Response2))
table(aux$Response2)
#Response3 vs MMR_sig
p<-ggplot(aux[aux$score_MMR_sig >0, ], 
          aes(x=Response3, y=score_MMR_sig)) + 
  geom_boxplot(outlier.shape = NA) + stat_compare_means() +
  annotate("text", label = paste('N = 27', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 4', sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_classic() +theme(legend.position="left")
gridExtra::grid.arrange(p, bottom="")
print(table(aux[aux$score_MMR_sig >0, ]$Response3))
table(aux$Response3)

#==============================================================================#

# Correlação MMR_sig vs TMB
p<-ggscatter(subset(aux, !(ID_Exoma %in% c("ROP-83", "ROP-107"))), x = "score_MMR_sig", y = "TMB",
             add = "reg.line",  # Add regressin line
             # add = "loess",  # Local Polynomial Regression Fitting
             #          cor.coef=T,
             cor.method = "spearman",
             add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
             conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",  label.x = .4)
gridExtra::grid.arrange(p, bottom="Spearman Test (67 samples)")

# Correlação MMR_sig vs MATH
p<-ggscatter(aux, x = "score_MMR_sig", y = "MATH_score",
             add = "reg.line",  # Add regressin line
             # add = "loess",  # Local Polynomial Regression Fitting
             #          cor.coef=T,
             cor.method = "spearman",
             add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
             conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",  label.x = .4,label.y = 70)
gridExtra::grid.arrange(p, bottom="Spearman Test (69 samples)")


#==============================================================================#
#dMMR vs MATH_score
p<-ggplot(sampleTableCount.Clinical, aes(x=dMMR, y=MATH_score)) + 
  geom_boxplot(outlier.shape = NA) + stat_compare_means() +
  annotate("text", label = paste('N = 57', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 12', sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_classic() +theme(legend.position="left")
gridExtra::grid.arrange(p, bottom="")
table(sampleTableCount.Clinical$dMMR)

#dMMR vs TMB
sampleTableCount.Clinical3 <- subset(sampleTableCount.Clinical, !(ID_Exoma %in% c("ROP-107", "ROP-83")))
p<-ggplot(sampleTableCount.Clinical3, aes(x=dMMR, y=TMB)) + 
  geom_boxplot(outlier.shape = NA) + stat_compare_means() +
  annotate("text", label = paste('N = 57', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 10', sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_classic() +theme(legend.position="left")
gridExtra::grid.arrange(p, bottom="(67 samples)")
table(sampleTableCount.Clinical3$dMMR)


#==============================================================================#
#                   DISTRIBUICAO DAS MAF NAS 69 AMOSTRAS  ####
#==============================================================================#
# Para plotar a distuibuiC'C#o dos MAFs estou usando as mutaC'C5es NON_SYNONIMOUS + SINONIMOUS
#somatic.69<-as.data.frame(fread(paste0(output,"/somatic_mutation.72samples.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))
# 
# aux1 <- select(subset(somatic.67, (SAMPLE %in% Clinical$Sample)  ), SAMPLE, MAF)
# aux2 <- as.data.frame(table(aux1), stringsAsFactors = F)
# aux2$MAF <- as.double(aux2$MAF)
# aux2 <- merge(aux2, Clinical[,c("Sample","ID_Exoma","Response1","Response2")], by.x= "SAMPLE", by.y = "Sample", all.x = T)

# MAFs das amostras outliers    ####
#==============================================================================#
# # nCRT-R Response1 
# P<- ggplot(data=subset(aux2, Response2 == "nCRT-R"),
#            aes(x=MAF, y=Freq, group = ID_Exoma, colour= Response2))+
#   geom_col(colour= "#08415c") +
#   # scale_y_sqrt()+
#   scale_x_continuous(breaks=seq(0,1,0.1))+
#   facet_wrap(vars(ID_Exoma))+theme_minimal()+
#   theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")
# P
# # "nCRT-NR Not-Metastatic"
# P<- ggplot(data=subset(aux2, Response2 == "nCRT-NR Not-Metastatic"),
#            aes(x=MAF, y=Freq, group = ID_Exoma, colour= Response2))+
#   geom_col(colour= "#ef767a") +
#   # scale_y_sqrt()+
#   scale_x_continuous(breaks=seq(0,1,0.1))+
#   facet_wrap(vars(ID_Exoma))+theme_minimal()+
#   theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")
# P
# # "nCRT-NR Metastatic"
# P<- ggplot(data=subset(aux2, Response2 == "nCRT-NR Metastatic"),
#            aes(x=MAF, y=Freq, group = ID_Exoma, colour= Response2))+
#   geom_col(colour= "#cc2956") +
#   # scale_y_sqrt()+
#   scale_x_continuous(breaks=seq(0,1,0.1))+
#   facet_wrap(vars(ID_Exoma))+theme_minimal()+
#   theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")
# P


#==============================================================================#
# Fisher Test  - Drivers genes, dMMR ####
#==============================================================================#

sampleTableCount.Clinical$DriverGene <-"notDriver"
mutados<- unique(subset(somatic.69, class == "driver", select = ID_Exoma))
sampleTableCount.Clinical[sampleTableCount.Clinical$ID_Exoma %in% mutados$ID_Exoma]$DriverGene <-"RectalDriver"

sampleTableCount.Clinical$dMMR <-"pMMR"
mutados<- unique(subset(somatic.69, class == "dMMR", select = ID_Exoma))
sampleTableCount.Clinical[sampleTableCount.Clinical$ID_Exoma %in% mutados$ID_Exoma]$dMMR <-"dMMR"

classeGenes <- c("DriverGene", "dMMR")
respostas <- c("Response1", "Response2", "Response3")

for (geneClass in classeGenes) {
  for (resp in respostas) {
    # resp <- "Response1"
    # geneClass<- "dMMR"
    dt3<- matrix(table(subset(sampleTableCount.Clinical, select = c(geneClass, resp))), nr=2)
    if(ncol(dt3)>1){
      resultado<-fisher.test(dt3)
      valor_p <- resultado$p.value
      if(valor_p <= 0.05){
        print("=============================================================")
        print(geneClass)
        print(resp)
        print(as.matrix(table(subset(sampleTableCount.Clinical, select = c(geneClass, resp))), nr=2))
        print(resultado)
        print("=============================================================")
        print(chisq.test(dt3))
        print("=============================================================")
      }
    }
  }
}


#==============================================================================#
# Fisher Test  - recurrentemente mutated genes ####
# number of patients with a specific gene mutated, associated with the outcomes
#==============================================================================#

# sink(file = paste0(output,"/Fisher_Test_dnds.mutGenes.FDR10.txt"), append = FALSE)

dndsGenes <-fread(paste0(output,"somatic.67.mut_signif_genes.tsv"), header = T, na.strings=c("NA"), fill=TRUE, check.names = FALSE)
# dndsGenes <- subset(dndsGenes, qglobal_cv <= 0.1, select = "gene_name")

for (resp in respostas) {
  for (gene in dndsGenes$gene_name) {
    # resp <- "Response1"
    # gene="FBXW7"
    aux <-sampleTableCount.Clinical2
    aux$mutGene <-"not_dNdSGene"
    mutados<- unique(subset(somatic.69,  Gene.refGene == gene, select = ID_Exoma))
    aux[aux$ID_Exoma %in% mutados$ID_Exoma, "mutGene"] <-"dNdSGene"
    aux1<- as.data.frame(unique(subset(aux, select = c("ID_Exoma", "mutGene", resp))))
    
    dt3<- matrix(table(subset(aux1, select = c("mutGene", resp))), nr=2)
    
    if(ncol(dt3)>1){
      resultado<-fisher.test(dt3)
      valor_p <- resultado$p.value
      if(valor_p <= 0.05){
        print("=============================================================")
        print(resp)
        print(gene)
        print(as.matrix(table(subset(aux1, select = c("mutGene", resp))), nr=2))
        print(resultado)
        print("=============================================================")
        print(chisq.test(dt3))
        print("=============================================================")
        print("=============================================================")
      }
    }
  }
}

# sink()

#==============================================================================#
# Fisher Test  - Drivers genes, dMMR ####
#==============================================================================#

# 
# #==============================================================================#
# #        Fisher e Chi.sq - Verificando a associação de alguma assinatura
# #   Teste para saber qual assinatura esta esta associada com o desfecho clinico
# #==============================================================================#
final_Ass <- fread("table_signatures.tsv")
aux<-as.data.frame(unique(subset(final_Ass,  sig_id==id, select=c("SAMPLE"))), stringsAsFactors = F)
List_Ass <-unique(final_Ass$sig_id)
# 
table_sig<-sampleTableCount.Clinical2
respostas <- c("Response1", "Response2", "Response3")

for (id in List_Ass) {
  for (resp in respostas) {
    # id="28"
    # resp= "Response1"
    aux<-as.data.frame(unique(subset(final_Ass,  sig_id==id, select=c("SAMPLE"))), stringsAsFactors = F)
    table_sig$MMR_sig = "0"
    table_sig[table_sig$SAMPLE %in% aux$SAMPLE]$MMR_sig <-"1"
    
    # dt3<- matrix(table(table_sig[,c("MMR_sig", "Response1")]), nr=2)
    dt3<- matrix(table(subset(table_sig, select = c("MMR_sig", resp))), nr=2)
    if(ncol(dt3)>1){
      resultado<-fisher.test(dt3)
      valor_p <- resultado$p.value
      if(valor_p <= 0.05){
        print("=============================================================")
        print(id)
        print(resp)
        print(as.matrix(table(subset(table_sig, select = c("MMR_sig", resp))), nr=2))
        print(resultado)
        print("=============================================================")
        print(chisq.test(dt3))
        print("=============================================================")
      }
    }
  }
}

# 
# #==============================================================================#
#          Chi.sq - Verificando a correlação de alguma assinatura
#           Considerando apenas as amostras com dMMR
#            SEM RESULTADOS SIGNIFICATIVOS
# #==============================================================================#
final_Ass <- fread("table_signatures.tsv")
aux<-as.data.frame(unique(subset(final_Ass,  sig_id==id, select=c("SAMPLE"))), stringsAsFactors = F)
List_Ass <-unique(final_Ass$sig_id)
# 
table_sig<-subset(sampleTableCount.Clinical, dMMR == "dMMR")
respostas <- c("Response1", "Response2", "Response3")

final_Ass<-subset(final_Ass, dMMR == "dMMR")
for (id in List_Ass) {
  for (resp in respostas) {
    # id="4"
    # resp= "Response2"
    aux<-as.data.frame(unique(subset(final_Ass,  sig_id==id, select=c("SAMPLE"))), stringsAsFactors = F)
    table_sig$MMR_sig = "0"
    table_sig[table_sig$SAMPLE %in% aux$SAMPLE]$MMR_sig <-"1"
    
    # dt3<- matrix(table(table_sig[,c("MMR_sig", "Response1")]), nr=2)
    dt3<- matrix(table(subset(table_sig, select = c("MMR_sig", resp))), nr=2)
    if(ncol(dt3)>1){
      resultado<-fisher.test(dt3)
      valor_p <- resultado$p.value
      if(valor_p <= 0.05){
        print("=============================================================")
        print(id)
        print(resp)
        print(as.matrix(table(subset(table_sig, select = c("MMR_sig", resp))), nr=2))
        print(resultado)
        print("=============================================================")
        print(chisq.test(dt3))
        print("=============================================================")
      }
    }
  }
}


#==============================================================================#
# CARREGANDO INFORMAÇÕES DE SVs
#==============================================================================#
tab_sv<- fread(file="../resultados_SV/tabela_count.SVs.tsv", header=T, sep="\t", fill=TRUE, check.names = FALSE)
tab_sv<-tab_sv %>% 
  mutate_all(replace_na, 0)
sampleTableCount.Clinical.SV <- sampleTableCount.Clinical

sampleTableCount.Clinical.SV <- merge(sampleTableCount.Clinical.SV, tab_sv, by.x = "ID_Exoma", by.y = "sample", all.x = T)
sampleTableCount.Clinical.SV$SV[is.na(sampleTableCount.Clinical.SV$SV)] <- 0
sampleTableCount.Clinical.SV$BND[is.na(sampleTableCount.Clinical.SV$BND)] <- 0
sampleTableCount.Clinical.SV$DEL[is.na(sampleTableCount.Clinical.SV$DEL)] <- 0
sampleTableCount.Clinical.SV$DUP[is.na(sampleTableCount.Clinical.SV$DUP)] <- 0
sampleTableCount.Clinical.SV$INV[is.na(sampleTableCount.Clinical.SV$INV)] <- 0
sampleTableCount.Clinical.SV$status_SV <- "sem SV"
sampleTableCount.Clinical.SV[sampleTableCount.Clinical.SV$SV > 0, "status_SV"] <- "com SV"


sampleTableCount.Clinical.SV$status_INV <- "sem_SV"
sampleTableCount.Clinical.SV$status_DUP <- "sem_SV"
sampleTableCount.Clinical.SV$status_DEL <- "sem_SV"
sampleTableCount.Clinical.SV$status_BND <- "sem_SV"

sampleTableCount.Clinical.SV[sampleTableCount.Clinical.SV$INV >0, ]$status_INV <- "com_SV"
sampleTableCount.Clinical.SV[sampleTableCount.Clinical.SV$DUP >0, ]$status_DUP <- "com_SV"
sampleTableCount.Clinical.SV[sampleTableCount.Clinical.SV$DEL >0, ]$status_DEL <- "com_SV"
sampleTableCount.Clinical.SV[sampleTableCount.Clinical.SV$BND >0, ]$status_BND <- "com_SV"

aux <-sampleTableCount.Clinical.SV
eventos_SV <- c("status_SV","status_INV","status_DUP","status_DEL","status_BND")

for (sv in eventos_SV) {
  for (resp in respostas) {
    # sv="status_SV"
    # resp= "Response1"
    
    dt3<- matrix(table(subset(sampleTableCount.Clinical.SV, select = c(sv, resp))), nr=2)
    dt3
    if(ncol(dt3)>1){
      resultado<-fisher.test(dt3)
      resultado1<-chisq.test(dt3)
      valor_p <- resultado$p.value
      valor_p1 <- resultado1$p.value
      if((valor_p <= 0.05)||(valor_p1 <= 0.05)){
        print("=============================================================")
        print(sv)
        print(resp)
        print(as.matrix(table(subset(aux, select = c(sv, resp))), nr=2))
        print("=============================================================")
        print(resultado)
        print("=============================================================")
        print(chisq.test(dt3))
        print("=============================================================")
      }
    }
  }
}






# =============================================================================#
#                        CURVA ROC  TMB - RESPOSTA
# =============================================================================#
library(pROC)

# dir.create("curva_roc")

# pdf("curva_roc/TMB_Resposta.pdf", width =7, height = 5)
#=========================================#
# resp.trat<-resp.trat[resp.trat$WXS == "Y",]
aux<-as.data.frame(sampleTableCount.Clinical2[,c('TMB','Response1')])
ggplot(aux, aes(x = TMB, fill = Response1)) +
  #geom_density(position = "identity", alpha = 0.5) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 20) +
  scale_fill_manual(values = c("blue", "red")) +
  ggtitle("TMB Distribuition (67 samples)") + theme_classic()+
  xlab("TMB") #+ ylab("Freq.")


# ============================#
# curva roc Completa vs Incompleta
# ============================#
aux<-as.data.frame(sampleTableCount.Clinical2[, c('TMB','Response1')],)
aux$Response1 <-ifelse(aux$Response1 == "nCRT-R", 1, 0) # USAR VITAL.STATUS OU RESPOSTA, RESPOSTA2

#=============================#
# Dividindo o dataset em treino e teste
# Criar o subset de treino
train <- aux %>% dplyr::sample_frac(.70)
# Criar o subset de teste com antijoin (pega tudo que não pertence)
test <- aux  

# Rodando o modelo
fit_reg_log <- 
  glm(Response1 ~ TMB,
      family = binomial(link = 'logit'),
      data = train)

# Fazendo as predições
pred_reg_log <- predict(fit_reg_log, newdata = test, type = "response")

# Organizando a tabela de dados para calcular as métricas da curva ROC
pred_roc_reg_log <- 
  dplyr::tibble(
    pred_reg_log,
    "survived" = as.factor(as.numeric(test$Response1)-1)
  ) %>% arrange(desc(pred_reg_log))

# Criando objeto com as métricas para curva ROC
roc_reg_log <- pROC::roc(pred_roc_reg_log$survived , pred_roc_reg_log$pred_reg_log, percent = TRUE)
# Se desejar, é possível (e bem simples) utilizar o próprio pacote pROC para plotar a curva ROC.
par(pty = "s")
plot.roc(
  main = "ROC Curve TMB (67 samples) \n nCRT-R vs nCRT-NR",
  roc_reg_log,
  print.auc = TRUE,
  legacy.axes = TRUE
  # xlab = "Taxa de Falso Positivo (100 - Especificidade)",
  # ylab = "Taxa de Verdadeiro Positivo (Sensibilidade)"
)



# =============================================================================#
#                        CURVA ROC  MATH_score - RESPOSTA
# =============================================================================#
library(pROC)
# 
# dir.create("curva_roc")
# 
# pdf("curva_roc/TMB_Resposta.pdf", width =7, height = 5)
#=========================================#
# resp.trat<-resp.trat[resp.trat$WXS == "Y",]

sampleTableCount.Clinical2.NotMetastatic = subset(sampleTableCount.Clinical2, Response2 != "nCRT-NR Metastatic")

aux<-as.data.frame(sampleTableCount.Clinical2.NotMetastatic[,c('MATH_score','Response1')])
ggplot(aux, aes(x = MATH_score, fill = Response1)) +
  #geom_density(position = "identity", alpha = 0.5) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 20) +
  scale_fill_manual(values = c("blue", "red")) +
  ggtitle("MATH_score Distribuition (67 samples)") + theme_classic()+
  xlab("MATH_score") #+ ylab("Freq.")

# ============================#
# curva roc Completa vs Incompleta
# ============================#
aux<-as.data.frame(sampleTableCount.Clinical2.NotMetastatic[, c('MATH_score','Response1')],)
aux$Response1 <-ifelse(aux$Response1 == "nCRT-R", 1, 0) # USAR VITAL.STATUS OU RESPOSTA, RESPOSTA2

#=============================#
# Dividindo o dataset em treino e teste
# Criar o subset de treino
train <- aux %>% dplyr::sample_frac(.70)
# Criar o subset de teste com antijoin (pega tudo que não pertence)
test <- aux  

# Rodando o modelo
fit_reg_log <- 
  glm(Response1 ~ MATH_score,
      family = binomial(link = 'logit'),
      data = train)

# Fazendo as predições
pred_reg_log <- predict(fit_reg_log, newdata = test, type = "response")

# Organizando a tabela de dados para calcular as métricas da curva ROC
pred_roc_reg_log <- 
  dplyr::tibble(
    pred_reg_log,
    "survived" = as.factor(as.numeric(test$Response1)-1)
  ) %>% arrange(desc(pred_reg_log))

# Criando objeto com as métricas para curva ROC
roc_reg_log <- pROC::roc(pred_roc_reg_log$survived , pred_roc_reg_log$pred_reg_log, percent = TRUE)
# Se desejar, é possível (e bem simples) utilizar o próprio pacote pROC para plotar a curva ROC.
par(pty = "s")
plot.roc(
  main = "ROC Curve MATH_score (51 samples)\n nCRT-R vs nCRT-NR",
  roc_reg_log,
  print.auc = TRUE,
  legacy.axes = TRUE
  # xlab = "Taxa de Falso Positivo (100 - Especificidade)",
  # ylab = "Taxa de Verdadeiro Positivo (Sensibilidade)"
)


