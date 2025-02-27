### VLIRA JUN2024 ####
### Step1 - Reestruturação dos dados de variantes identificados pelo Mutect2 GATK-v4.6 ###
# ==============================================================================#
# Nessa versão do script eu seleciono apenas as 72 amostras
# 


# Na versão anterior (2023-09-26) faltava excluir mutações na lista de HLA
# Apliquei os filtros:
# gnomAD: 0.005
# Abraom: 0.005
# shared >=2 & cosmic <=3
# VAF < 0.10
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

# Fisher Test (dnds_genes, MMR, drivers)
# Calcula Curva ROC do TMB e MATH_score para Resposta1

#==============================================================================#

library(data.table)
library(stringr)

library(dplyr)
library(reshape2)
library(tidyr)

library(ggplot2)
library(devtools)
library(dndscv)
library(ggpubr)
library(xlsx)
library(maftools)
library(pheatmap)
library(pROC)


wd <- "D:/Github-Projects/Pipeline_Gatk-Mutect2/Result_Mutect2.GATK4.6.2014-06-14/"
setwd(wd)
dir.create("output-Step1")
output <- "output-Step1/"

#==============================================================================#
# Dataset reorganization ####
#==============================================================================#
all_variants <- fread(file="Final_GATK.4.6_annotated.txt", header=T, sep="\t", na.strings=c(".","NA"), fill=TRUE, check.names = FALSE)

all_variants<-all_variants[all_variants$Chr != "Chr"]
all_variants <- all_variants[, -c("Otherinfo1","Otherinfo2","Otherinfo3")]

# Rename some columns from table for clarity
colnames(all_variants)<- c("SAMPLE",colnames(all_variants[,2:37]),"CHR","POS","ID","REF","ALT","QUAL","FILTER", "INFO", "FORMAT","GENOTYPE")

# Creating INDEX (CHR,POS,REF,ALT) and select the main columns
index <- paste(all_variants$CHR, all_variants$POS, all_variants$REF, all_variants$ALT, sep="_")
all_variants <- cbind(index, all_variants)

# Format the SAMPLE info 
all_variants <- separate(data= all_variants, col=SAMPLE, into=c("SAMPLE"), sep = ".hg38_multianno.txt") 
all_variants$EXOME_ID <- all_variants$SAMPLE
all_variants <- separate(data= all_variants, col=EXOME_ID, into=c("EXOME_ID"), sep = "-ExC85")
nrow(table(all_variants$SAMPLE))
#==============================================================================#
fwrite(all_variants, paste0(output,"/all_variants.allSamples.tsv"), quote = F, sep="\t")

## Remover amostras clinicas excluidas ###
samples.excl <- c("ROP-116","ROP-118","ROP-120","ROP-121","ROP-122","ROP-123","ROP-125","ROP-126","ROP-127","ROP-130","ROP-131","ROP-132","ROP-23","ROP-84")
all_variants.72sample <- subset(all_variants, (!EXOME_ID %in% samples.excl))
length(unique(all_variants.72sample$EXOME_ID)) # 72
length(unique(all_variants.72sample$index))   # 79558

## Separate genotype info ("GT","AD","AF","DP","F1R2","F2R1","FAD") ###
unique(all_variants.72sample$FORMAT) #"GT:AD:AF:DP:F1R2:F2R1:FAD:SB"    "GT:AD:AF:DP:F1R2:F2R1:FAD:PGT:PID:PS:SB"
FORMAT<-strsplit(unique(all_variants.72sample$FORMAT), ":")[[2]]
FORMAT

# Separate genotype information
all_variants.72sample<- separate(data = all_variants.72sample, col = GENOTYPE, into=FORMAT, sep = ":")

## Separate Alleles ###
all_variants.72sample$aux = all_variants.72sample$GT
all_variants.72sample$aux = str_replace(all_variants.72sample$aux,"/","|")   
all_variants.72sample<-cbind(all_variants.72sample, read.table(text = as.character(all_variants.72sample$aux), sep = "|"))
names(all_variants.72sample)[c(ncol(all_variants.72sample)-1,ncol(all_variants.72sample))] <- c('GT1','GT2')
all_variants.72sample<- subset(all_variants.72sample, select= -c(aux))
setDT(all_variants.72sample)[,paste0("AD", 1:2) := tstrsplit(AD, ",")]
table(all_variants.72sample$GT)
# 0/1   0|1   1|0 
# 95475 18530  1572 

## CALCULAR VAF e COV ###
all_variants.72sample$COV <- as.numeric(all_variants.72sample$AD1)+as.numeric(all_variants.72sample$AD2)
all_variants.72sample$VAF <- round(as.numeric(all_variants.72sample$AD2)/(as.numeric(all_variants.72sample$AD1)+as.numeric(all_variants.72sample$AD2)),5)

# Quantificando mutações compartilhadas >=2
aux <- data.frame(table(all_variants.72sample$index))
all_variants.72sample <- merge(all_variants.72sample, aux, by.x = "index", by.y = "Var1")
colnames(all_variants.72sample)[names(all_variants.72sample)=="Freq"]<-"MUTATION_shared"
nrow(unique(all_variants.72sample[,"index"]))                   #79558
#==============================================================================#
all_var.72<- all_variants.72sample

## Quantify mutations present on COSMIC 
#==============================================================================#
all_var.72 <- separate(data = all_var.72, col = cosmic99, into = c("COSMIC_ID","OCCURENCE"), sep = "OCCURENCE=")
all_var.72$COSMIC_OCCURENCE <- all_var.72$OCCURENCE
all_var.72 <- separate(data =all_var.72, col =OCCURENCE, into = c("V_1","V_2","V_3","V_4","V_5"), sep = ",")
aux <- as.data.frame(lapply(all_var.72[,c("V_1","V_2","V_3","V_4","V_5")], function(y) as.numeric(gsub("\\(.*\\)", "", y))))
aux[is.na(aux)] <- 0
all_var.72 <- cbind(subset(all_var.72, select = -c(V_1,V_2,V_3,V_4,V_5)), aux)
all_var.72$COSMIC_count <- rowSums(all_var.72[,c("V_1","V_2","V_3","V_4","V_5")])
all_var.72 <- dplyr::select(all_var.72, !starts_with("V_"))

## Classify mutations like (none, driver, dMMMR) 
#==============================================================================#
driver<-fread("../input_Somatic-Filter/driver-genes.txt", header = F, col.names = "gene")
MMR<-fread("../input_Somatic-Filter/MMR_genes.txt", header = F, col.names = "gene")
# HR <-fread("../input_Somatic-Filter/MMR_genes.txt", header = F, col.names = "gene")
# BER<-fread("../input_Somatic-Filter/MMR_genes.txt", header = F, col.names = "gene")
# NER<-fread("../input_Somatic-Filter/MMR_genes.txt", header = F, col.names = "gene")
# NHEJ <-fread("../input_Somatic-Filter/MMR_genes.txt", header = F, col.names = "gene")
# dnaRepair<-fread("../input_Somatic-Filter/MMR_genes.txt", header = F, col.names = "gene")
all_var.72$class = "none"
all_var.72$class[all_var.72$Gene.refGene %in% driver$gene] <-"driver"
all_var.72$class[all_var.72$Gene.refGene %in% MMR$gene] <-"MMR"

all_var.72$driver = 0
all_var.72$MMR = 0
all_var.72$POLE = 0
all_var.72$driver[all_var.72$Gene.refGene %in% driver$gene] <- 1
all_var.72$MMR[all_var.72$Gene.refGene %in% MMR$gene] <- 1
all_var.72$POLE[all_var.72$Gene.refGene == "POLE"] <- 1

fwrite(all_var.72, paste0(output,"/all_variants.72sample.tsv"), quote = F, sep="\t")

#==============================================================================#
#                    STEP 4 - APPLY SOMATIC FILTER        ####
#==============================================================================#
all<-fread(paste0(output,"all_variants.72sample.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)

length(unique(all$index)) # 79558
## Quantifica MUTAÇOES GERMINATIVAS 
Germline_gnomad <- which(all$gnomad41_exome_AF_grpmax > 0.0005)
nrow(unique(all[Germline_gnomad,"index"]))  # 5385; 8723 nao-redundate germinativos GnomAD
## Quantifica MUTAÇOES GERMINATICAS ABRAOM 
Germline_abraom <- which(all$abraom_freq > 0.0005)
nrow(unique(all[Germline_abraom,"index"])) # 3848; 7006 nao-redundate germinativos Abraom

## Quantifica MUTAÇOES not-rare gnomAD
not_rare2exclude <- which(all$gnomad41_exome_AF_grpmax > 0.0005)
length(not_rare2exclude)                      # 10165; 15193 excluidos
nrow(unique(all[not_rare2exclude,"index"]))   # 5385; 8723 nao-redundate excluidos
mut.GF1 <- all[-c(not_rare2exclude),]
length(unique(mut.GF1$index))                 # 74173; 70835 variantes ficaram

## Quantifica MUTAÇOES not-rare no Abraom
not_rareAbraom2exclude <- which(mut.GF1$abraom_freq > 0.0005)
length(not_rareAbraom2exclude)                            # 387; 2127 excluidos
nrow(unique(mut.GF1[not_rareAbraom2exclude,"index"]))     # 184;1811 nao-redundate excluidos
mut.GF2 <- mut.GF1[-c(not_rareAbraom2exclude),]
length(unique(mut.GF2$index))                             #73989; 69024 variantes restantes
somatic_Var <- mut.GF2

# Selecionando ads Germinativa
germline_Var <- setdiff(all, somatic_Var)

length(unique(germline_Var$index))              # 5569; 10534 variantes germinativas 
fwrite(germline_Var, paste0(output,"/germline_variants.tsv"), quote = F, sep="\t")

length(unique(somatic_Var$index))               # 73989; 69024 variantes somaticas
fwrite(somatic_Var, paste0(output,"/somatic_variants.tsv"), quote = F, sep="\t")
#==============================================================================#
#CHECK WHITELIST
whitelist <- fread(file= paste0("D:/PROJETOS-HSL_BP/hg38lft_genome_whitelist.bed"))
colnames(whitelist) <- c("Chr","Start","End", "Gene")

tabela1 <- germline_Var[,c("Chr","Start","End","index","VAF","COV", "AAChange.refGene","gnomad41_exome_AF_grpmax","abraom_freq")]

is_in_whitelist <- function(value, chr, whitelist) {
  any(whitelist$Chr == chr & whitelist$Start <= value & value <= whitelist$End)
}
table1_checked <- tabela1 %>%
  rowwise() %>%
  mutate(
    in_range_col1 = is_in_whitelist(Start, Chr, whitelist),
    in_range_col2 = is_in_whitelist(End, Chr, whitelist)
  ) %>%
  ungroup()

# NÃO EXISTE VARIANTES GERMINATIVAS QUE ESTÃO EM REGIÕES DA WHITELIST
View(subset(table1_checked, (in_range_col1 ==T) | (in_range_col2 == T)))

## FILTRO: REMOVE mutações compartilhadas >=2 e <= 3 ocorencias no COSMIC ####
#==============================================================================#
somatic_Var<-as.data.frame(fread(paste0(output,"somatic_variants.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))
length(unique(somatic_Var$index)) # 73989; 69024 variantes somaticas

shared2exclude <- which( somatic_Var$MUTATION_shared >= 2 & somatic_Var$COSMIC_count <= 3)
length(unique(somatic_Var[shared2exclude,"index"]))   # 9011; 8456
somatic_Var.Filt_ShareCosmic <- somatic_Var[-shared2exclude,]  
length(unique(somatic_Var.Filt_ShareCosmic[,"index"]))          #64978;60568    
length(unique(somatic_Var.Filt_ShareCosmic[,"SAMPLE"]))          #72    

somatic_Var.Filt_ShareCosmic$COSMIC_count[somatic_Var.Filt_ShareCosmic$COSMIC_count>= 5] <- paste(somatic_Var.Filt_ShareCosmic$COSMIC_count[somatic_Var.Filt_ShareCosmic$COSMIC_count>= 5], "+")

fwrite(somatic_Var.Filt_ShareCosmic, paste0(output,"/somatic_variants.FiltSharedCosmic.tsv"), quote = F, sep="\t")
somatic_Var.Filt_SC<-as.data.frame(fread(paste0(output,"somatic_variants.FiltSharedCosmic.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))

#==============================================================================#
#                   DISTRIBUICAO DAS VAF NAS 72 AMOSTRAS  ####
#==============================================================================#
Clinical <- xlsx::read.xlsx("D:/PROJETOS-HSL_BP/Dados Clinicos 2024 ALLversions WXS-TCRseq.xlsx", sheetName = "73samples", header = T)
Clinical <- subset(Clinical, select=-Paciente)
Clinical <- Clinical[Clinical$WXS== "Y",]
Clinical <- Clinical[rowSums(is.na(Clinical)) != ncol(Clinical),]

somatic_Var.Filt_SC<-as.data.frame(fread(paste0(output,"somatic_variants.FiltSharedCosmic.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))
length(unique(somatic_Var.Filt_SC$index)) # 64978; 60568
aux1 <- select(subset(somatic_Var.Filt_SC, (SAMPLE %in% Clinical$Sample)  ), SAMPLE, VAF)
aux2 <- as.data.frame(table(aux1), stringsAsFactors = F)
aux2$VAF <- as.double(aux2$VAF)
aux2 <- merge(aux2, Clinical[,c("Sample","EXOME_ID","Response1","Response2")], by.x= "SAMPLE", by.y = "Sample", all.x = T)

## FILTRO: REMOVE MUTACOES EM HLA ####
#==============================================================================#
genesHLA <- fread("../input_Somatic-Filter/genesHLA.txt", col.names = "gene", header = F) 
all.hla<- subset(somatic_Var.Filt_SC, (Gene.refGene %in% genesHLA$gene))
length(unique(all.hla$index))    # 55; 52

hla2exclude <- which(somatic_Var.Filt_SC$Gene.refGene %in% genesHLA$gene) 
length(unique(somatic_Var.Filt_SC[hla2exclude,"index"]))  #  52
somatic_Var.Filt_SC_HLA <- somatic_Var.Filt_SC[-hla2exclude,]
length(unique(somatic_Var.Filt_SC_HLA$index))  #64923; 60516

## FILTRO: REMOVE MUTACOES  MINOR ALLELIC FREQ. < 0.10 ####
#==============================================================================#
VAF2exclude <- which(somatic_Var.Filt_SC_HLA$VAF<0.1) 
length(unique(somatic_Var.Filt_SC_HLA[VAF2exclude,"index"]))  #  26404; 25224 excluidos
somatic_Var.Filt_SC_HLA_VAF <- somatic_Var.Filt_SC_HLA[-VAF2exclude,]
length(unique(somatic_Var.Filt_SC_HLA_VAF$index))  #38589; 35346

## FILTRO: REMOVE MUTACOES COV < 30 ####
#==============================================================================#
cov2exclude <- which(somatic_Var.Filt_SC_HLA_VAF$COV<30) 
length(unique(somatic_Var.Filt_SC_HLA_VAF[cov2exclude,"index"]))  #  1504 excluidas
somatic_Var.Filt_SC_HLA_VAF_COV <- somatic_Var.Filt_SC_HLA_VAF[-cov2exclude,]
length(unique(somatic_Var.Filt_SC_HLA_VAF_COV$index))  # 33849 restante


fwrite(somatic_Var.Filt_SC_HLA_VAF_COV, paste0(output,"/somatic_variants.Filters.tsv"), quote = F, sep="\t")

somatic_Var.Filters<-as.data.frame(fread(paste0(output,"somatic_variants.Filters.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))
## FILTRO: SELECT EXONIC AND SPLICING ####
#==============================================================================#
somatic_Var.CodRegion <- subset(somatic_Var.Filters, Func.refGene %in% c("splicing", "exonic", "exonic;splicing", "splicing;exonic"))
length(unique(somatic_Var.CodRegion$index))  #37571;34395
unique(somatic_Var.CodRegion$Func.refGene)
#34126765 - xGen Exome Research Panel V2 Target Regions
Mb <- 34126765/1000000    
#Mb=34.126
tmb <- as.data.frame(table(somatic_Var.CodRegion$SAMPLE))
names(tmb) <- c("SAMPLE", "N_Mutations_TMB")
tmb$TMB <- round(tmb$N_Mutations/Mb,2)

fwrite(tmb, paste0(output,"/TMB.tsv"), quote = F, sep="\t")
fwrite(somatic_Var.CodRegion, paste0(output,"/somatic_variants.CodRegion.tsv"), quote = F, sep="\t")


#==============================================================================#
somatic_Var.CodRegion<-as.data.frame(fread(paste0(output,"somatic_variants.CodRegion.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))
length(unique(somatic_Var.CodRegion[,"index"]))     # 33693

aux<-(unique(somatic_Var.CodRegion[,c("index", "ExonicFunc.refGene")]))  
as.data.frame(table(aux$ExonicFunc.refGene))

#  REMOVE MUTACOES CLASSIFICADAS COM SINONIMAS (synonymous SNV)   8134
#==============================================================================#
somatic = subset(somatic_Var.CodRegion, !(ExonicFunc.refGene %in% c("synonymous SNV")))
length(unique(somatic[,"index"]))     # 25559 restantes

### SALVAR A TABELA COM MUTACOES SOMATICAS  ####
#==============================================================================#
fwrite(somatic, paste0(output,"/somatic_mutation.tsv"), quote = F, sep="\t")
#==============================================================================#
somatic<-as.data.frame(fread(paste0(output,"/somatic_mutation.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))

#==============================================================================#
# CALCULAR MATH SCORE  ####
#==============================================================================#
Tab.toMATH <- somatic
somatic.LENGTH <- as.data.frame(tapply(Tab.toMATH$VAF, Tab.toMATH$SAMPLE, length))
somatic.LENGTH$SAMPLE <- row.names(somatic.LENGTH)
names(somatic.LENGTH) <- c("N_mutations", "SAMPLE")
somatic.MEDIAN <- as.data.frame(tapply(Tab.toMATH$VAF, Tab.toMATH$SAMPLE, median))
somatic.MEDIAN$SAMPLE <-row.names(somatic.MEDIAN)
names(somatic.MEDIAN) <- c("VAF_median", "SAMPLE")
somatic.MAD <- as.data.frame(tapply(Tab.toMATH$VAF, Tab.toMATH$SAMPLE, mad))
somatic.MAD$SAMPLE <-row.names(somatic.MAD)
names(somatic.MAD) <- c("VAF_MAD", "SAMPLE")
match_somatic <- match(somatic.MAD$SAMPLE, somatic.MEDIAN$SAMPLE)
somatic.MAD$VAF_median <- somatic.MEDIAN$VAF_median[match_somatic]
match_somatic_length <- match(somatic.MAD$SAMPLE, somatic.LENGTH$SAMPLE)
somatic.MAD$N_mutations <- somatic.LENGTH$N_mutations[match_somatic_length]
somatic.MATH <- somatic.MAD
somatic.MATH$MATH_score <- 100*(somatic.MATH$VAF_MAD/somatic.MATH$VAF_median)

#==============================================================================#
# TABELA COM SCORES DE HETEROGENIEDADE  ####
#==============================================================================#
somatic_Samples<- merge(somatic.MATH, tmb, by = "SAMPLE")

# SALVAR A TABELA MATH-SCORE, TMB DE CADA AMOSTRA  ###
#==============================================================================#
fwrite(somatic_Samples, paste0(output,"/MATH_TMB-SAMPLES.tsv"), quote = F, sep="\t")

# SALVAR ARQUIVOS BED PARA FILTRAR VCF ####
#==============================================================================#
# samples<-unique(Clinical$Sample)
# dir.create(paste0(output,"mutations_to_filter_fromVCF/"))
# for (sp in samples) {
#   BED<-subset(somatic, SAMPLE==sp, select = c("CHR", "POS","REF","ALT","index"))
#   BED<-unique(BED)
#    fwrite(BED,paste0(output,"mutations_to_filter_fromVCF/",sp,".bed"), sep = "\t", col.names = T, row.names = F, quote = F)
# }

# SALVAR ARQUIVOS BED COM SINONIMAS PARA FILTRAR VCF ####
#==============================================================================#
samples<-unique(Clinical$Sample)
dir.create(paste0(output,"mutations_to_filter_fromVCF/"))
for (sp in samples) {
  BED<-subset(somatic_Var.CodRegion, SAMPLE==sp, select = c("CHR", "POS","REF","ALT","index"))
  BED<-unique(BED)
  fwrite(BED,paste0(output,"mutations_to_filter_fromVCF/",sp,".synonymous.bed"), sep = "\t", col.names = T, row.names = F, quote = F)
}

## Carregando dados de MATH_SCORE, TMB, MUT, SAMPLE
#sampleTableCount<- fread(paste0(output,"/MATH_TMB-SAMPLES.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)
sampleTableCount<-somatic_Samples

aux<-subset(somatic_Samples, SAMPLE %in% Clinical$Sample)
fwrite(aux, paste0(output,"/MATH_TMB-72.SAMPLES.tsv"), quote = F, sep="\t")

## Filtrando mutações dos  pacientes finais 
sampleTableCount<-subset(sampleTableCount, (SAMPLE %in% Clinical$Sample) )

## Juntando tabelas TMB_MATH_#MUT +  dados_clinicos 
sampleTableCount.Clinical <- merge(sampleTableCount, Clinical, by.x = "SAMPLE", by.y = "Sample")
sampleTableCount.Clinical2 <- subset(sampleTableCount.Clinical, !(EXOME_ID %in% c("ROP-107", "ROP-83", "ROP-90")), )

fwrite(sampleTableCount.Clinical, paste0(output,"/sampleTableCount.Clinical.tsv"), quote = F, sep="\t")
fwrite(sampleTableCount.Clinical2, paste0(output,"/sampleTableCount.Clinical2.tsv"), quote = F, sep="\t")






#==============================================================================#
# PIPELINE EXECUTADO CORRETAMENTE ATE AQUI,
# POIS EU NAO RODEI A ANALISE DE ASSINATURA MUTACIONAL (MUTALISK)
#==============================================================================#







#==============================================================================#
#                 MMR vs MMR_sig vs MATH_score vs TMB                      ####
#==============================================================================#

mutsig <-fread("table_signatures.tsv", header = T)
mutsig <- merge(Clinical[,c(1,2)], mutsig, by.x= "Sample",by.y = "SAMPLE")

aux<-unique(subset(somatic, class == "MMR",select = c("SAMPLE","class") ))
sampleTableCount.Clinical$MMR = "pMMR"
sampleTableCount.Clinical$MMR[sampleTableCount.Clinical$SAMPLE %in% aux$SAMPLE] <-'MMR'

aux <- as.data.frame(unique(subset(mutsig, cosmic_sig == "MMR_sig", select = "Sample")))
sampleTableCount.Clinical$MMR_Sig = "others_sig"
sampleTableCount.Clinical[sampleTableCount.Clinical$SAMPLE %in% aux$Sample, "MMR_Sig"] <- "MMR_sig"

aux <- sampleTableCount.Clinical

print(table(aux$MMR))
print(table(aux$MMR_Sig))

print(as.matrix(table(aux[,c("MMR", "MMR_Sig")]), nr=2))
dt3<- matrix(table(aux[,c("MMR", "MMR_Sig")]), nr=2)
print(fisher.test(dt3))
print(chisq.test(dt3))

score_MMR_sig <- subset(mutsig, cosmic_sig=="MMR_sig") %>% group_by(Sample, Response1, Response3) %>% 
  dplyr::summarise(score_MMR_sig = sum(perc),
                   .groups = 'drop')

# aux$score_MMR_sig <- 0
aux<- merge(aux, score_MMR_sig[,c(1,4)], by.x= "SAMPLE", by.y= "Sample", all.x = T)
aux$score_MMR_sig[is.na(aux$score_MMR_sig)] <- 0

p<-ggplot(subset(aux, score_MMR_sig>0), aes(x=MMR, y=score_MMR_sig)) + 
  geom_boxplot(aes(colour = MMR),outlier.shape = NA) + stat_compare_means() +
  annotate("text", label = paste('N = 25', sep = ''), x = 2, y = -0.05, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 6', sep = ''), x = 1, y = -0.05, angle = 0, size = 3.75) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_classic() +theme(legend.position="none")
gridExtra::grid.arrange(p, bottom="")
print(table(aux[aux$score_MMR_sig >0, ]$MMR))
print(table(aux$MMR))

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
p<-ggscatter(subset(aux, !(EXOME_ID %in% c("ROP-83", "ROP-107"))), x = "score_MMR_sig", y = "TMB",
             add = "reg.line",  # Add regressin line
             # add = "loess",  # Local Polynomial Regression Fitting
             #          cor.coef=T,
             cor.method = "spearman",
             add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
             conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",  label.x = .4)
gridExtra::grid.arrange(p, bottom="Spearman Test (70 samples)")

# Correlação MMR_sig vs MATH
p<-ggscatter(aux, x = "score_MMR_sig", y = "MATH_score",
             add = "reg.line",  # Add regressin line
             # add = "loess",  # Local Polynomial Regression Fitting
             #          cor.coef=T,
             cor.method = "spearman",
             add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
             conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",  label.x = .4,label.y = 70)
gridExtra::grid.arrange(p, bottom="Spearman Test (72 samples)")


#==============================================================================#
#MMR vs MATH_score
p<-ggplot(sampleTableCount.Clinical, aes(x=MMR, y=MATH_score)) + 
  geom_boxplot(outlier.shape = NA) + stat_compare_means() +
  annotate("text", label = paste('N = 57', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 12', sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_classic() +theme(legend.position="left")
gridExtra::grid.arrange(p, bottom="")
table(sampleTableCount.Clinical$MMR)

#MMR vs TMB
sampleTableCount.Clinical3 <- subset(sampleTableCount.Clinical, !(EXOME_ID %in% c("ROP-107", "ROP-83")))
p<-ggplot(sampleTableCount.Clinical3, aes(x=MMR, y=TMB)) + 
  geom_boxplot(outlier.shape = NA) + stat_compare_means() +
  annotate("text", label = paste('N = 57', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 10', sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_classic() +theme(legend.position="left")
gridExtra::grid.arrange(p, bottom="(70 samples)")
table(sampleTableCount.Clinical3$MMR)



#==============================================================================#
# Fisher Test / Chisq Test - Drivers genes, MMR ####
#==============================================================================#

sampleTableCount.Clinical$DriverGene <-"notDriver"
mutados<- unique(subset(somatic.72, class == "driver", select = EXOME_ID))
sampleTableCount.Clinical[sampleTableCount.Clinical$EXOME_ID %in% mutados$EXOME_ID,]$DriverGene <-"RectalDriver"

sampleTableCount.Clinical$MMR <-"pMMR"
mutados<- unique(subset(somatic.72, class == "MMR", select = EXOME_ID))
sampleTableCount.Clinical[sampleTableCount.Clinical$EXOME_ID %in% mutados$EXOME_ID,]$MMR <-"MMR"

classeGenes <- c("DriverGene", "MMR")
respostas <- c("Response1", "Response2", "Response3")

for (geneClass in classeGenes) {
  for (resp in respostas) {
    dt3<- matrix(table(subset(sampleTableCount.Clinical, select = c(geneClass, resp))), nr=2)
    if(ncol(dt3)>1){
      resultado<-fisher.test(dt3)
      resultado1<-chisq.test(dt3)
      valor_p <- resultado$p.value
      valor_p1 <- resultado1$p.value
      if((valor_p <= 0.05)||(valor_p1 <= 0.05)){
        print("=============================================================")
        print("=============================================================")
        print(geneClass)
        print(resp)
        print(as.matrix(table(subset(sampleTableCount.Clinical, select = c(geneClass, resp))), nr=2))
        print(resultado)
        print("========== Chisq Test =======================================")
        print(chisq.test(dt3))
        print("========== Fisher Test ======================================")
        print(fisher.test(dt3))
      }
    }
  }
}


#==============================================================================#
# Fisher Test  - recurrentemente mutated genes ####
# number of patients with a specific gene mutated, associated with the outcomes
#==============================================================================#

sink(file = paste0(output,"/Fisher_Test_dnds.mutGenes.FDR10.txt"), append = FALSE)
dndsGenes <-fread(paste0(output,"somatic.72.mut_signif_genes.tsv"), header = T, na.strings=c("NA"), fill=TRUE, check.names = FALSE)
dndsGenes <- subset(dndsGenes, qglobal_cv <= 0.1, select = "gene_name")
for (resp in respostas) {
  for (gene in dndsGenes$gene_name) {
    aux <-sampleTableCount.Clinical
    aux$mutGene <-"not_dNdSGene"
    mutados<- unique(subset(somatic.72,  Gene.refGene == gene, select = EXOME_ID))
    aux[aux$EXOME_ID %in% mutados$EXOME_ID, "mutGene"] <-"dNdSGene"
    aux1<- as.data.frame(unique(subset(aux, select = c("EXOME_ID", "mutGene", resp))))
    dt3<- matrix(table(subset(aux1, select = c("mutGene", resp))), nr=2)
    if(ncol(dt3)>1){
      resultado<-fisher.test(dt3)
      resultado1<-chisq.test(dt3)
      valor_p <- resultado$p.value
      valor_p1 <- resultado1$p.value
      if((valor_p <= 0.05)||(valor_p1 <= 0.05)){
        print("=============================================================")
        print("=============================================================")
        print(resp)
        print(gene)
        print(as.matrix(table(subset(aux1, select = c("mutGene", resp))), nr=2))
        print(resultado)
        print("========== Chisq Test =======================================")
        print(chisq.test(dt3))
        print("========== Fisher Test ======================================")
        print(fisher.test(dt3))
      }
    }
  }
}

sink()

#==============================================================================#
# Fisher Test  - Drivers genes, MMR ####
#==============================================================================#
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
#           Considerando apenas as amostras com MMR
#            SEM RESULTADOS SIGNIFICATIVOS
# #==============================================================================#
final_Ass <- fread("table_signatures.tsv")
aux<-as.data.frame(unique(subset(final_Ass,  sig_id==id, select=c("SAMPLE"))), stringsAsFactors = F)
List_Ass <-unique(final_Ass$sig_id)
# 
table_sig<-subset(sampleTableCount.Clinical, MMR == "MMR")
respostas <- c("Response1", "Response2", "Response3")

final_Ass<-subset(final_Ass, MMR == "MMR")
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

sampleTableCount.Clinical.SV <- merge(sampleTableCount.Clinical.SV, tab_sv, by.x = "EXOME_ID", by.y = "sample", all.x = T)
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
  ggtitle("TMB Distribuition (70 samples)") + theme_classic()+
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
  main = "ROC Curve TMB (70 samples) \n nCRT-R vs nCRT-NR",
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
  ggtitle("MATH_score Distribuition (70 samples)") + theme_classic()+
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



### codigo utilizado para plotar figuras utilizadas no labmeeting de 18/06/2024



tab <- read.xlsx(paste0(output, "/Table.72.summary.xlsx"), sheet = "Clin.SNV.CNV.SV", )
p <- ggplot(tab, aes(x=Response2, y=GAIN, fill=Response2)) + 
  geom_boxplot(outlier.shape = NA)+
  stat_compare_means()+
  scale_y_sqrt()+
  scale_fill_manual(values = c("#08415c","#cc2956", "#ef767a"), breaks = c("nCRT-R", "nCRT-NR Metastatic", "nCRT-NR Not-Metastatic")) +
  annotate("text", label = paste('N = 16', sep = ''), x = 3, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 38', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 18', sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  stat_compare_means(comparisons = my_comparisons)+
  theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="CNAs (GAIN) correlation Response2 (72 samples)")
table(tab$Response2)

p <- ggplot(tab, aes(x=Response2, y=LOSS, fill=Response2)) + 
  geom_boxplot(outlier.shape = NA)+
  stat_compare_means()+
  scale_y_sqrt()+
  scale_fill_manual(values = c("#08415c","#cc2956", "#ef767a"), breaks = c("nCRT-R", "nCRT-NR Metastatic", "nCRT-NR Not-Metastatic")) +
  annotate("text", label = paste('N = 16', sep = ''), x = 3, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 38', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 18', sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  stat_compare_means(comparisons = my_comparisons)+
  theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="CNAs (LOSS) correlation Response2 (72 samples)")
table(tab$Response2)

p <- ggplot(tab, aes(x=Response2, y=Mean.AF, fill=Response2)) + 
  geom_boxplot(outlier.shape = NA)+
  stat_compare_means()+
  scale_fill_manual(values = c("#08415c","#cc2956", "#ef767a"), breaks = c("nCRT-R", "nCRT-NR Metastatic", "nCRT-NR Not-Metastatic")) +
  annotate("text", label = paste('N = 16', sep = ''), x = 3, y = 0.15, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 38', sep = ''), x = 2, y = 0.15, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 18', sep = ''), x = 1, y = 0.15, angle = 0, size = 3.75) +
  stat_compare_means(comparisons = my_comparisons)+
  theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="Mean.AF correlation Response2 (72 samples)")
table(tab$Response2)


#= SVs

p <- ggplot(tab, aes(x=Response2, y=MantaBND, fill=Response2)) + 
  geom_boxplot(outlier.shape = NA)+
  stat_compare_means()+
  scale_fill_manual(values = c("#08415c","#cc2956", "#ef767a"), breaks = c("nCRT-R", "nCRT-NR Metastatic", "nCRT-NR Not-Metastatic")) +
  annotate("text", label = paste('N = 16', sep = ''), x = 3, y = -1.5, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 38', sep = ''), x = 2, y = -1.5, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 18', sep = ''), x = 1, y = -1.5, angle = 0, size = 3.75) +
  stat_compare_means(comparisons = my_comparisons)+
  theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="MantaBND correlation Response2 (72 samples)")
table(tab$Response2)

p <- ggplot(tab, aes(x=Response2, y=MantaDEL, fill=Response2)) + 
  geom_boxplot(outlier.shape = NA)+
  stat_compare_means()+
  scale_fill_manual(values = c("#08415c","#cc2956", "#ef767a"), breaks = c("nCRT-R", "nCRT-NR Metastatic", "nCRT-NR Not-Metastatic")) +
  annotate("text", label = paste('N = 16', sep = ''), x = 3, y = 0.15, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 38', sep = ''), x = 2, y = 0.15, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 18', sep = ''), x = 1, y = 0.15, angle = 0, size = 3.75) +
  stat_compare_means(comparisons = my_comparisons)+
  theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="MantaDEL correlation Response2 (72 samples)")
table(tab$Response2)

p <- ggplot(tab, aes(x=Response2, y=MantaDUP, fill=Response2)) + 
  geom_boxplot(outlier.shape = NA)+
  stat_compare_means()+
  scale_fill_manual(values = c("#08415c","#cc2956", "#ef767a"), breaks = c("nCRT-R", "nCRT-NR Metastatic", "nCRT-NR Not-Metastatic")) +
  annotate("text", label = paste('N = 16', sep = ''), x = 3, y = 0.15, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 38', sep = ''), x = 2, y = 0.15, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 18', sep = ''), x = 1, y = 0.15, angle = 0, size = 3.75) +
  stat_compare_means(comparisons = my_comparisons)+
  theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="MantaDUP correlation Response2 (72 samples)")
table(tab$Response2)

p <- ggplot(tab, aes(x=Response2, y=MantaINS, fill=Response2)) + 
  geom_boxplot(outlier.shape = NA)+
  stat_compare_means()+
  scale_fill_manual(values = c("#08415c","#cc2956", "#ef767a"), breaks = c("nCRT-R", "nCRT-NR Metastatic", "nCRT-NR Not-Metastatic")) +
  annotate("text", label = paste('N = 16', sep = ''), x = 3, y = 0.15, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 38', sep = ''), x = 2, y = 0.15, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 18', sep = ''), x = 1, y = 0.15, angle = 0, size = 3.75) +
  stat_compare_means(comparisons = my_comparisons)+
  theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="MantaINS correlation Response2 (72 samples)")
table(tab$Response2)



