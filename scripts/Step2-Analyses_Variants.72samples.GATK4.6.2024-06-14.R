
library(data.table)
library(stringr)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(ggpubr)




wd <- "D:/Github-Projects/Pipeline_Gatk-Mutect2/Result_Mutect2.GATK4.6.2014-06-14/"
setwd(wd)
dir.create("output-Step1")
output <- "output-Step1/"
dir.create(paste0(output,"Figures"))

#==============================================================================#
somatic<-as.data.frame(fread(("output-Step1/somatic_mutation.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))
somatic_Samples<- fread(paste0("output-Step1/MATH_TMB-SAMPLES.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)
somatic_Var.CodRegion<-as.data.frame(fread(paste0("output-Step1/somatic_variants.CodRegion.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))

#==============================================================================#
# Carregando tabela de dados Clinicos ####
#==============================================================================#
Clinical <- xlsx::read.xlsx("D:/PROJETOS-HSL_BP/Dados Clinicos 2024 ALLversions WXS-TCRseq.xlsx", sheetName = "73samples", header = T)
Clinical <- subset(Clinical, select=-Paciente)
Clinical <- Clinical[Clinical$WXS== "Y",]
Clinical <- Clinical[rowSums(is.na(Clinical)) != ncol(Clinical),]


# SALVAR ARQUIVOS BED COM SINONIMAS PARA FILTRAR VCF ####
#==============================================================================#
samples<-unique(Clinical$Sample)
dir.create(paste0(output,"mutations_to_filter_fromVCF/"))
for (sp in samples) {
  BED<-subset(somatic_Var.CodRegion, SAMPLE==sp, select = c("CHR", "POS","REF","ALT","index"))
  BED<-unique(BED)
  fwrite(BED,paste0(output,"mutations_to_filter_fromVCF/",sp,".synonymous.bed"), sep = "\t", col.names = T, row.names = F, quote = F)
}

sampleTableCount<-somatic_Samples

## Juntando tabelas TMB_MATH_#MUT +  dados_clinicos 
sampleTableCount.Clinical <- merge(sampleTableCount, Clinical, by.x = "SAMPLE", by.y = "Sample")
sampleTableCount.Clinical2 <- subset(sampleTableCount.Clinical, !(EXOME_ID %in% c("ROP-107", "ROP-83", "ROP-90")), )


#==============================================================================#
### PLOT MATH_score ####
#==============================================================================#
somatic_Var.Filt_SC<-as.data.frame(fread(paste0(output,"somatic_variants.FiltSharedCosmic.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))

aux1<-as.data.frame(table(somatic_Var.Filt_SC$COV), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)

pdf(file = paste0(output,"Figures/Cov.somatic_Var.Filt_SC.pdf"), width = 9, height = 5)
p<- ggplot(data=head(aux1, 500), aes(x=Var1, y=Freq))+
  geom_col() + 
  geom_vline(xintercept= 20, linetype="dashed", color = "red") +
  xlab("Deapth of Coverage") +theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
LEGEND<- c("LEGEND: The figure display the Coverage distribution range < 500x, from somatic variants in 72 samples.
Only shared&Cosmic filter applied. The redline indicates threashold at 20x.")
gridExtra::grid.arrange(p, bottom=LEGEND)
dev.off()

aux1<-as.data.frame(table(somatic_Var.Filt_SC$VAF), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
pdf(file = paste0(output,"Figures/VAF.somatic_Var.Filt_SC.pdf"), width = 9, height = 5)
p<-ggplot(data=aux1, aes(x=Var1, y=Freq))+
  geom_col(width = .005) + 
  geom_vline(xintercept= 0.10, linetype="dashed", color = "red")+
  xlab("VAF") + theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
  LEGEND<- c("LEGEND: The figure display the VAF distribution from somatic variants in 72 samples.
  Only shared&Cosmic filter applied. The redline indicates threashold at VAF(< 0.1) ")
gridExtra::grid.arrange(p, bottom=LEGEND)
dev.off()

## MATH-SCORE colour by Resposta1 
pdf(file = paste0(output,"Figures/MATH-score.distribution.72samples.pdf"), width = 12, height = 6)
p<-ggplot(data=sampleTableCount.Clinical, aes(x= reorder(EXOME_ID, -MATH_score), y=MATH_score, fill=Response1)) +
  geom_bar(stat="identity")+theme_classic() +
  scale_fill_manual(values = c("#08415c", "#cc2936"), breaks = c("nCRT-R", "nCRT-NR")) +
  geom_hline(yintercept= median(sampleTableCount.Clinical$MATH_score), linetype="dashed", color = "black")+
  annotate(geom="text", x=50, y=50, label=paste("median MATH= ",sprintf("%.2f",median(sampleTableCount.Clinical$MATH_score))), color="black")+
  xlab("SAMPLES") +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gridExtra::grid.arrange(p, bottom="MATH-score distribution from 72 samples")
dev.off()


# Generate a bar plot displaying distinct ExonicFunc mutation types
aux<-as.data.frame(unique(somatic_Var.CodRegion[,c("index", "ExonicFunc.refGene")]))
aux <- as.data.frame(table(aux$ExonicFunc.refGene))
aux <- aux %>%
  arrange(desc(Freq)) %>%
  mutate(Var1 = factor(Var1, levels = unique(Var1)))
colnames(aux) <- c("ExonicFunc", "count")
# Create a horizontal bar plot with values in front of bars
# Supondo que seu dataframe se chama 'aux'
aux <- aux %>%  mutate(proportion = count / sum(count) * 100)
# Criar o gráfico com contagem e proporção

pdf(file = paste0(output,"Figures/ExonicFunc.refGene.distribution.72samples.pdf"), width = 8, height = 6)
p<-ggplot(data = aux, aes(x = ExonicFunc, y = count, fill = ExonicFunc)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = paste(count, sprintf("(%.1f%%)", proportion)), y = 18000),
            hjust = 0.1, color = "black") +
  coord_flip() +  theme_minimal() +
  labs(title = "Count of Distinct ExonicFunc Mutation Types",
       x = "ExonicFunc Mutation Type",
       y = "Count") +
  scale_fill_grey(start = 0.8, end = 0.2) +
  theme(legend.position = "none")
gridExtra::grid.arrange(p, bottom="ExonicFunc.refGene distribution from 72 samples")
dev.off()

#==============================================================================#
# PLOT BOXPLOT MATH_score CORRELATION CLINICAL
#==============================================================================#
# RESPONSE 1
table(sampleTableCount.Clinical$Response1)
pdf(file = paste0(output,"Figures/BOXPLOTS_MATH_score.72samples.pdf"), width = 8, height = 6)
p <- ggplot(sampleTableCount.Clinical  , aes(x=Response1, y=MATH_score, fill=Response1)) + 
  geom_boxplot(outlier.shape = NA)+ 
  scale_fill_manual(values = c("#08415c", "#cc2936"), breaks = c("nCRT-R", "nCRT-NR")) +
  annotate("text", label = paste('N = 16', sep = ''), x = 2, y = 3, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 56', sep = ''), x = 1, y = 3, angle = 0, size = 3.75) +
  stat_compare_means()+ theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
table(sampleTableCount.Clinical$Response1)
gridExtra::grid.arrange(p, bottom="MATH-score correlation nCRT-R x nCRT-NR (72 samples)")
#dev.off()

# sem Metastaticos
#==============================================================================#
aux<-sampleTableCount.Clinical[sampleTableCount.Clinical$Response3 != "Metastatic", ] 
table(aux[aux$Response2 != "nCRT-NR Metastatic", "Response2"])
#pdf(file = paste0(output,"Figures/boxplot2_MATH_score.72samples.pdf"), width = 8, height = 6)
p <- ggplot( aux, aes(x=Response1, y=MATH_score, fill=Response1)) + 
  geom_boxplot(outlier.shape = NA)+ 
  scale_fill_manual(values = c("darkgrey", "lightgrey"), breaks = c("nCRT-R", "nCRT-NR")) +
  annotate("text", label = paste('N = 15', sep = ''), x = 2, y = 3, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 38', sep = ''), x = 1, y = 3, angle = 0, size = 3.75) +
  stat_compare_means()+ theme_classic()+
  theme(legend.position="none") + geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="MATH-score correlation nCRT-R x nCRT-NR (not Metastatic)")
table(aux$Response1)
#dev.off()

my_comparisons <- list( c("nCRT-R", "nCRT-NR Not-Metastatic"), 
                        c("nCRT-NR Not-Metastatic", "nCRT-NR Metastatic"), 
                        c("nCRT-R", "nCRT-NR Metastatic") )
# RESPONSE 2
#pdf(file = paste0(output,"Figures/boxplot3_MATH_score.72samples.pdf"), width = 8, height = 6)
p <- ggplot(sampleTableCount.Clinical, aes(x=Response2, y=MATH_score, fill=Response2)) + 
  geom_boxplot(outlier.shape = NA)+
  stat_compare_means(label.y = 80)+
  scale_fill_manual(values = c("#08415c","#cc2956", "#ef767a"), breaks = c("nCRT-R", "nCRT-NR Metastatic", "nCRT-NR Not-Metastatic")) +
  annotate("text", label = paste('N = 16', sep = ''), x = 3, y = 3, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 38', sep = ''), x = 2, y = 3, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 18', sep = ''), x = 1, y = 3, angle = 0, size = 3.75) +
  stat_compare_means(comparisons = my_comparisons)+
  theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="MATH-score correlation Response2 (72 samples)")
table(sampleTableCount.Clinical$Response2)
dev.off()


## BARPLOT TMB ####
#==============================================================================#
## TMB colour by Resposta1 
pdf(file = paste0(output,"Figures/TMB_DISTRIBUTION.72samples.pdf"), width = 12, height = 6)
p<- ggplot(data=sampleTableCount.Clinical, aes(x= reorder(EXOME_ID, -TMB), y=TMB, fill=Response1)) +
  geom_bar(stat="identity")+theme_classic() +
  scale_y_sqrt()+
  scale_fill_manual(values = c("#08415c", "#cc2936"), breaks = c("nCRT-R", "nCRT-NR")) +
  geom_hline(yintercept= median(sampleTableCount.Clinical$TMB), linetype="dashed", color = "black")+
  annotate(geom="text", x=50, y=50, label=paste("median TMB= ",sprintf("%.2f",median(sampleTableCount.Clinical$TMB))), color="black")+
  xlab("SAMPLES")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gridExtra::grid.arrange(p, bottom="TMB distribution from 72 samples")

#==============================================================================#
# PLOT BOXPLOT TMB CORRELATION CLINICAL
#==============================================================================#
# RESPONSE 1
#pdf(file = paste0(output,"Figures/BOXPLOT1_TMB.72samples.pdf"), width = 12, height = 6)
p <- ggplot(sampleTableCount.Clinical  , aes(x=Response1, y=TMB, fill=Response1)) + 
  geom_boxplot(outlier.shape = NA)+ 
  scale_y_sqrt()+
  scale_fill_manual(values = c("#08415c", "#cc2936"), breaks = c("nCRT-R", "nCRT-NR")) +
  annotate("text", label = paste('N = 16', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 56', sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  stat_compare_means()+ theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="TMB correlation Response1 (72 samples)")
#dev.off()
table(sampleTableCount.Clinical$Response1)

# sem ROP-107, ROP-83 & ROP-83
#==============================================================================#
#pdf(file = paste0(output,"Figures/BOXPLOT1_TMB.69samples.pdf"), width = 12, height = 6)
p <- ggplot(sampleTableCount.Clinical2, aes(x=Response1, y=TMB, fill=Response1)) + 
  geom_boxplot(outlier.shape = NA)+ 
  # scale_y_sqrt()+
  scale_fill_manual(values = c("#08415c", "#cc2936"), breaks = c("nCRT-R", "nCRT-NR")) +
  annotate("text", label = paste('N = 13', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 56', sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  stat_compare_means()+ theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="TMB correlation Response1 (69 samples)")
#dev.off()
table(sampleTableCount.Clinical2$Response1)

my_comparisons <- list( c("nCRT-R", "nCRT-NR Not-Metastatic"), 
                        c("nCRT-NR Not-Metastatic", "nCRT-NR Metastatic"), 
                        c("nCRT-R", "nCRT-NR Metastatic") )
# RESPONSE 2
sampleTableCount.Clinical2<-sampleTableCount.Clinical[!(sampleTableCount.Clinical$EXOME_ID %in% c("ROP-107","ROP-83","ROP-90")),]

#pdf(file = paste0(output,"Figures/BOXPLOT2_TMB.72samples.pdf"), width = 12, height = 6)
p <- ggplot(sampleTableCount.Clinical, aes(x=Response2, y=TMB, fill=Response2)) + 
  geom_boxplot(outlier.shape = NA)+
  scale_y_sqrt()+
  stat_compare_means(label.y = 10)+
  stat_compare_means(comparisons = my_comparisons)+ 
  scale_fill_manual(values = c("#08415c","#cc2956", "#ef767a"), breaks = c("nCRT-R", "nCRT-NR Metastatic", "nCRT-NR Not-Metastatic")) +
  annotate("text", label = paste('N = 16', sep = ''), x = 3, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 38', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 18', sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="TMB correlation Response2 (72 samples)")
#dev.off()
table(sampleTableCount.Clinical$Response2)
# sem ROP-107 e ROP-83
#==============================================================================#
#pdf(file = paste0(output,"Figures/BOXPLOT2_TMB.69samples.pdf"), width = 12, height = 6)
p <- ggplot(sampleTableCount.Clinical2, aes(x=Response2, y=TMB, fill=Response2)) + 
  geom_boxplot(outlier.shape = NA)+
  scale_y_sqrt()+
  stat_compare_means(label.y = 5)+
  stat_compare_means(comparisons = my_comparisons)+ 
  scale_fill_manual(values = c("#08415c","#cc2956", "#ef767a"), breaks = c("nCRT-R", "nCRT-NR Metastatic", "nCRT-NR Not-Metastatic")) +
  annotate("text", label = paste('N = 13', sep = ''), x = 3, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 38', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 18', sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="TMB correlation Response2 (69 samples)")
#dev.off()
table(sampleTableCount.Clinical2$Response2)

# RESPONSE 3
#pdf(file = paste0(output,"Figures/BOXPLOT3_TMB.72samples.pdf"), width = 12, height = 6)
p <- ggplot(sampleTableCount.Clinical, aes(x=Response3, y=TMB, fill=Response3)) + 
  geom_boxplot(outlier.shape = NA)+
  scale_y_sqrt()+
  stat_compare_means(label.y = 10)+
  scale_fill_manual(values = c("#cc2956", "#ef764a"), breaks = c("Metastatic", "Not-Metastatic")) +
  annotate("text", label = paste('N = 53', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 19', sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="TMB correlation Response3 (72 samples)")
#dev.off()
table(sampleTableCount.Clinical$Response3)
# sem ROP-107 e ROP-83
#==============================================================================#
#pdf(file = paste0(output,"Figures/BOXPLOT3_TMB.69samples.pdf"), width = 12, height = 6)
p <- ggplot(sampleTableCount.Clinical2, aes(x=Response3, y=TMB, fill=Response3)) + 
  geom_boxplot(outlier.shape = NA)+
  # scale_y_sqrt()+
  stat_compare_means(label.y = 10)+
  scale_fill_manual(values = c("#cc2956", "#ef764a"), breaks = c("Metastatic", "Not-Metastatic")) +
  annotate("text", label = paste('N = 50', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 19', sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  theme_classic()+theme(legend.position="none")+
  geom_jitter(shape=16, position=position_jitter(0.2))
gridExtra::grid.arrange(p, bottom="TMB correlation Response3 (69 samples)")
dev.off()

table(sampleTableCount.Clinical2$Response3)



#==============================================================================#
#                 CORRELACAO TMB vs MATH-SCORE           ####
#==============================================================================#
pdf(file = paste0(output,"Figures/tmb_x_math.pdf"), width = 12, height = 6)

shapiro.test(sampleTableCount.Clinical$TMB)
qqnorm(sampleTableCount.Clinical$TMB)
qqline(sampleTableCount.Clinical$TMB)

shapiro.test(sampleTableCount.Clinical$MATH_score)
qqnorm(sampleTableCount.Clinical$MATH_score)
qqline(sampleTableCount.Clinical$MATH_score)

sampleTableCount.Clinical2 <- subset(sampleTableCount.Clinical, !(EXOME_ID %in% c("ROP-107", "ROP-83", "ROP-90")), )

shapiro.test(sampleTableCount.Clinical2$TMB)
qqnorm(sampleTableCount.Clinical2$TMB)
qqline(sampleTableCount.Clinical2$TMB)

shapiro.test(sampleTableCount.Clinical2$MATH_score)
qqnorm(sampleTableCount.Clinical2$MATH_score)
qqline(sampleTableCount.Clinical2$MATH_score)

p<-ggscatter(sampleTableCount.Clinical, x = "MATH_score", y = "TMB",
             #add = "reg.line",  # Add regressin line
              add = "loess",  # Local Polynomial Regression Fitting
             #          cor.coef=T,
             cor.method = "spearman",
             add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
             conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",  label.x = 30,label.y = 100)
gridExtra::grid.arrange(p, bottom="Spearman Test (72 samples)")

sp <- ggscatter(sampleTableCount.Clinical2, x = "MATH_score", y = "TMB",
                #add = "reg.line",  # Add regressin line
                 add = "loess",  # Local Polynomial Regression Fitting
                # cor.coef=T,
                cor.method = "spearman",
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 30, label.y = 17)
gridExtra::grid.arrange(sp, bottom="Spearman Test (69 samples)")

#sem Metastatic
sp <- ggscatter(subset(sampleTableCount.Clinical2,Response3 == "Not-Metastatic"), x = "MATH_score", y = "TMB",
                #add = "reg.line",  # Add regressin line
                add = "loess",  # Local Polynomial Regression Fitting
                cor.method = "spearman",
                color = "white",
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 40)+ 
  geom_point(aes(color = Response2))+
  theme(legend.position="top")
gridExtra::grid.arrange(sp, bottom="Spearman Test (50 samples) excluded Metastatic")
nrow(subset(sampleTableCount.Clinical2,Response3 == "Not-Metastatic", ))

# RESPONSE 2  excluded metastatic
sp <- ggscatter(subset(sampleTableCount.Clinical2,Response2 != "nCRT-NR Metastatic"), x = "MATH_score", y = "TMB",
                add = "reg.line",  # Add regressin line
                cor.method = "spearman",
                color = "white",
                add.params = list(color = "Response2", fill = "Response2"),
               # add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = F # Add confidence interval
) + stat_cor(method = "spearman",  aes(color = Response2), label.x = 40)+
  geom_point(aes(color = Response2))+
  theme(legend.position="top")
gridExtra::grid.arrange(sp, bottom="Spearman Test (50 samples) excluded Metastatic")

# RESPONSE 1
sp <- ggscatter(sampleTableCount.Clinical2, x = "MATH_score", y = "TMB",
                #add = "reg.line",  # Add regressin line
                add = "loess",  # Local Polynomial Regression Fitting
                cor.method = "spearman",
                color = "white",
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman",  aes(color = Response1), label.x = 40)+
  geom_point(aes(color = Response1))+
  theme(legend.position="top")
gridExtra::grid.arrange(sp, bottom="Spearman Test (69 samples)")
nrow(sampleTableCount.Clinical2)

# RESPONSE 2
sp <- ggscatter(sampleTableCount.Clinical2, x = "MATH_score", y = "TMB",
                #add = "reg.line",  # Add regressin line
                add = "loess",  # Local Polynomial Regression Fitting
                cor.method = "spearman",
                color = "white",
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman",  aes(color = Response2), label.x = 40)+
  geom_point(aes(color = Response2))+
  theme(legend.position="top")
gridExtra::grid.arrange(sp, bottom="Spearman Test (69 samples)")

# RESPONSE 2 NOVO
sp <- ggscatter(sampleTableCount.Clinical2, x = "MATH_score", y = "TMB",
                add = "reg.line",  # Add regressin line
                #add = "loess",  # Local Polynomial Regression Fitting
                cor.method = "spearman",
                color = "white",
                add.params = list(color = "Response2", fill = "Response2"), # Customize reg. line
                conf.int = F # Add confidence interval
) + stat_cor(method = "spearman",  aes(color = Response2,fill = "Response2"), label.x = 40)+
  geom_point(aes(color = Response2))+
  theme(legend.position="top")
gridExtra::grid.arrange(sp, bottom="Spearman Test (69 samples)")

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
gridExtra::grid.arrange(sp, bottom="Spearman Test (70 samples)")

dev.off()


#==============================================================================#
# TOP 50 GENES MAIS MUTADOS ####
#==============================================================================#
driver<-fread("../input_Somatic-Filter/driver-genes.txt", header = F, col.names = "gene")
MMR<-fread("../input_Somatic-Filter/MMR_genes.txt", header = F, col.names = "gene")
mut.genes=unique(dplyr::select(somatic, c("index","Gene.refGene")))
aux2 <- as.data.frame(table(mut.genes$Gene.refGene))
names(aux2) <- c("Gene", "Freq_Mutation")
aux2 <- (aux2[order(aux2$Freq, decreasing = T),])
aux2$class <-'none'
aux2$class[aux2$Gene %in% driver$gene] <-'driver'
aux2$class[aux2$Gene %in% MMR$gene] <-'MMR'
#fwrite(aux2, "output_Somatic-Filter/countMut_by_Genes.tsv", quote = F, sep="\t")
aux2=head(aux2,100)
top50<- aux2$Gene

pdf(file = paste0(output,"Figures/countMut_Genes.pdf"), width = 12, height = 6)
p<-ggplot(data=aux2, aes(x= reorder(Gene, -Freq_Mutation), y=Freq_Mutation, fill=class)) +
  geom_bar(stat="identity" )+
  # scale_y_sqrt()+
  theme_classic() +xlab("Gene")
p +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



