library(reshape2)
library(maftools)
library(data.table)
library(ggpubr)


wd <- "D:/Github-Projects/Pipeline_Gatk-Mutect2/scripts/subset_notMetastatic/"
setwd(wd)
input<- "D:/Github-Projects/Pipeline_Gatk-Mutect2/Result_Mutect2.GATK4.6.2014-06-14/"
dir.create("output-SomaticFilter")
output <- "output-SomaticFilter/"
dir.create(paste0(output,"Figures"))

#==============================================================================#
# GERAR TABELA NO FORMATO DO Maftools  ####
#==============================================================================#
Clinical <- xlsx::read.xlsx("D:/PROJETOS-HSL_BP/Dados Clinicos 2024 ALLversions WXS-TCRseq.xlsx", sheetName = "73samples", header = T)
Clinical <- subset(Clinical, select=-Paciente)
Clinical <- Clinical[Clinical$WXS== "Y",]
Clinical <- Clinical[rowSums(is.na(Clinical)) != ncol(Clinical),]
# seleciona amostras não metastaticas
Clinical <- Clinical[Clinical$Response3 == "Not-Metastatic",]
Clinical$Tumor_Sample_Barcode <- Clinical$EXOME_ID

#==============================================================================#
# GERAR TABELA NO FORMATO DO Maftools  ####
#==============================================================================#
syn.Mutation.MAF<-as.data.frame(fread(paste0(output,"syn.Mutation.MAF.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))
nonsynList<-unique(syn.Mutation.MAF[syn.Mutation.MAF$Variant_Classification != "Silent",]$Variant_Classification)
maf.Tab <- read.maf(maf = syn.Mutation.MAF
                    , # vc_nonSyn = nonsynList, 
                    clinicalData = Clinical)

length(unique(syn.Mutation.MAF$index))
length(unique(syn.Mutation.MAF$index))

#==============================================================================#
# MATH score - ITH  ####
#==============================================================================#
rop.math <- math.score(maf = maf.Tab, vafCol = 'VAF' )
# SALVAR A TABELA MATH-SCORE  ###
#==============================================================================#
fwrite(rop.math, paste0(output,"/rop.math.tsv"), quote = F, sep="\t")

#==============================================================================#
# CONSTRUINDO SAMPLES_DATA ####
#==============================================================================#
Samples_Data.3<- merge(Samples_Data.2, rop.math, by= "Tumor_Sample_Barcode" )

# rop.het = inferHeterogeneity(maf = maf.Tab, vafCol = 'VAF')
# # #Visualizing results. Highlighting those variants on copynumber altered variants.
# plotClusters(clusters = rop.het, genes=dnds_genes , showCNvars = TRUE)

#==============================================================================#
# COMPARA  MATH SCORE  Maftool vs Implementação MANUAL ####
#==============================================================================#
nonsyn_MUT <- maf.Tab@data
#nonsyn_MUT <- syn.Mutation.MAF[!(syn.Mutation.MAF$Variant_Classification %in% c("Silent")),]

somatic.LENGTH <- as.data.frame(tapply(nonsyn_MUT$VAF, nonsyn_MUT$Tumor_Sample_Barcode, length))
somatic.LENGTH$Tumor_Sample_Barcode <- row.names(somatic.LENGTH)
names(somatic.LENGTH) <- c("N_mutations", "Tumor_Sample_Barcode")
somatic.MEDIAN <- as.data.frame(tapply(nonsyn_MUT$VAF, nonsyn_MUT$Tumor_Sample_Barcode, median))
somatic.MEDIAN$Tumor_Sample_Barcode <-row.names(somatic.MEDIAN)
names(somatic.MEDIAN) <- c("VAF_median", "Tumor_Sample_Barcode")
somatic.MAD <- as.data.frame(tapply(nonsyn_MUT$VAF, nonsyn_MUT$Tumor_Sample_Barcode, mad))
somatic.MAD$Tumor_Sample_Barcode <-row.names(somatic.MAD)
names(somatic.MAD) <- c("VAF_MAD", "Tumor_Sample_Barcode")
match_somatic <- match(somatic.MAD$Tumor_Sample_Barcode, somatic.MEDIAN$Tumor_Sample_Barcode)
somatic.MAD$VAF_median <- somatic.MEDIAN$VAF_median[match_somatic]
match_somatic_length <- match(somatic.MAD$Tumor_Sample_Barcode, somatic.LENGTH$Tumor_Sample_Barcode)
somatic.MAD$N_mutations <- somatic.LENGTH$N_mutations[match_somatic_length]
somatic.MATH <- somatic.MAD
somatic.MATH$MATH_score <- 100*(somatic.MATH$VAF_MAD/somatic.MATH$VAF_median)

MATHs_Compare<- base::merge(somatic.MATH, rop.math, by = "Tumor_Sample_Barcode")#
MATHs_Compare

#==============================================================================#
# COV
#==============================================================================#
aux1<-as.data.frame(table(syn.Mutation.MAF$COV), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)

pdf(file = paste0(output,"Figures/Cov.syn.Mutation.pdf"), width = 9, height = 5)
p<- ggplot(data=head(aux1, 500), aes(x=Var1, y=Freq))+
  geom_col() + 
  geom_vline(xintercept= 30, linetype="dashed", color = "red") +
  xlab("Deapth of Coverage") +theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
LEGEND<- c("LEGEND: The figure display the Coverage distribution range < 500x, from somatic variants in 53 samples.
 The redline indicates threashold at 30x.")
gridExtra::grid.arrange(p, bottom=LEGEND)
dev.off()

#VAF
#==============================================================================#
aux1<-as.data.frame(table(syn.Mutation.MAF$VAF), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
pdf(file = paste0(output,"Figures/VAF.syn.Mutation.MAF.pdf"), width = 9, height = 5)
p<-ggplot(data=aux1, aes(x=Var1, y=Freq))+
  geom_col(width = .005) + 
  geom_vline(xintercept= 0.10, linetype="dashed", color = "red")+
  xlab("VAF") + theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
LEGEND<- c("LEGEND: The figure display the VAF distribution from somatic variants in 53 samples.
  The redline indicates threashold at VAF(< 0.1) ")
gridExtra::grid.arrange(p, bottom=LEGEND)
dev.off()


#==============================================================================#
### PLOT MATH_score ####
#==============================================================================#

paleta <- c("darkgrey", "lightgrey")
paleta <- c("#08415c", "#cc2936")

# Verificar normalidade
shapiro.test(Samples_Data$MATH[dados$Response1 == "nCRT-R"])
shapiro.test(Samples_Data$MATH[dados$Response1 == "nCRT-NR"])
# Verificar homogeneidade das variâncias
leveneTest(MATH ~ Response1, data = Samples_Data)

## MATH-SCORE colour by Resposta1 
pdf(file = paste0(output,"Figures/MATH_Score.Plots.pdf"), width = 10, height = 8)
p<-ggplot(data=Samples_Data, aes(x= reorder(Tumor_Sample_Barcode, -MATH), y=MATH, fill=Response1)) +
  geom_bar(stat="identity")+theme_classic() + theme(legend.position="top")+
  scale_fill_manual(values = paleta, breaks = c("nCRT-R", "nCRT-NR")) +
  geom_hline(yintercept= median(Samples_Data$MATH), linetype="dashed", color = "black")+
  annotate(geom="text", x=40, y=60, label=paste("median MATH= ",sprintf("%.2f",median(Samples_Data$MATH))), color="black")+
  xlab("SAMPLES") +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gridExtra::grid.arrange(p, bottom="MATH-score distribution")
#dev.off()
#pdf(file = paste0(output,"Figures/boxplot2_MATH_BOXPLOT.pdf"), width = 8, height = 6)
count_N <- table(Samples_Data$Response1)
p <- ggplot( Samples_Data, aes(x=Response1, y=MATH, fill=Response1)) + 
  geom_boxplot(outlier.shape = NA)+ 
  scale_fill_manual(values = paleta, breaks = c("nCRT-R", "nCRT-NR")) +
  annotate("text", label = paste0('N = ',count_N[[2]]), x = 2, y = 3, angle = 0, size = 3.75) +
  annotate("text", label = paste0('N = ',count_N[[1]]), x = 1, y = 3, angle = 0, size = 3.75) +
  stat_compare_means()+ theme_classic()+
  theme(legend.position="none") + geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="MATH-score correlation nCRT-nR x nCRT-R ")
table(Samples_Data$Response1)
dev.off()



