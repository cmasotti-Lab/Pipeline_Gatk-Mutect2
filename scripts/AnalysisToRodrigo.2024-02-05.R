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


color.R<-"gray"
color.NR<-"white"

wd <- "D:/PROJETOS-HSL_BP/resultados_Mutect2-2024/"
# wd <- "/media/vande/HDD/PROJETOS-HSL_BP/resultados_Mutect2-2023/"
setwd(wd)
dir.create("output-Mutect2_SomaticFilters.2024-02-05/")
dir.create("output-Mutect2_SomaticFilters.2024-02-05/plots")
output <- "output-Mutect2_SomaticFilters.2024-02-05/"

Clinical69<-fread(paste0(output,"/Clinical.69.tsv"), quote = F, sep="\t")
Clinical69$Resposta3 <- ifelse(Clinical69$ID_Exoma == 'ROP-50', 'metastatico', Clinical69$Resposta3)
Clinical69$Response3 <- ifelse(Clinical69$ID_Exoma == 'ROP-50', 'Metastatic', Clinical69$Response3)
#fwrite(Clinical69,paste0(output,"/Clinical.69.tsv"), quote = F, sep="\t")

somatic.69 <- fread(paste0(output,"/somatic_mutation.69samples.tsv"), quote = F, sep="\t")

sampleTableCount.Clinical <- fread(paste0(output,"/sampleTableCount.Clinical.69.tsv"), quote = F, sep="\t")
sampleTableCount.Clinical$Resposta3 <- ifelse(sampleTableCount.Clinical$ID_Exoma == 'ROP-50', 'metastatico', sampleTableCount.Clinical$Resposta3)
sampleTableCount.Clinical$Response3 <- ifelse(sampleTableCount.Clinical$ID_Exoma == 'ROP-50', 'Metastatic', sampleTableCount.Clinical$Response3)
#fwrite(sampleTableCount.Clinical, paste0(output,"/sampleTableCount.Clinical.69.tsv"), quote = F, sep="\t" )

sampleTableCount.Clinical$ClinicalOutcome <- sampleTableCount.Clinical$Response1
# EXCLUIDO OS METASTATICOS
sampleTableCount.Clinical <- subset(sampleTableCount.Clinical, Response3 !=  "Metastatic")
sampleTableCount.Clinical2 <- subset(sampleTableCount.Clinical, !(ID_Exoma %in% c("ROP-107", "ROP-83")), )

somatic.52 <- subset(somatic.69, SAMPLE %in% sampleTableCount.Clinical$SAMPLE )

immuno<-fread(paste0(output,"/ROP_immuno_69pat_final.tsv"), quote = F, sep="\t")
immuno<-subset(immuno, Sample_Name %in% sampleTableCount.Clinical$ID_Exoma, 
               select = c('Sample_Name','Clones','Heterozygosis','HED_Score','Mean_class','NAL','nal_bin','NAL_snvs','NAL_indel','HLA.mut'))

sampleTableCount.Clinical$DriverGene <-"notDriver"
mutados<- unique(subset(somatic.52, class == "driver", select = ID_Exoma))
sampleTableCount.Clinical[sampleTableCount.Clinical$ID_Exoma %in% mutados$ID_Exoma]$DriverGene <-"RectalDriver"

sampleTableCount.Clinical$dMMR <-"pMMR"
mutados<- unique(subset(somatic.52, class == "dMMR", select = ID_Exoma))
sampleTableCount.Clinical[sampleTableCount.Clinical$ID_Exoma %in% mutados$ID_Exoma]$dMMR <-"dMMR"

#==============================================================================#
### PLOT MATH_score ####
#==============================================================================#
table(sampleTableCount.Clinical$ClinicalOutcome)
aux<-(table(sampleTableCount.Clinical$ClinicalOutcome))
p<-ggplot(data=sampleTableCount.Clinical, aes(x= reorder(ID_Exoma, -MATH_score), y=MATH_score, fill=ClinicalOutcome)) +
  geom_bar(stat="identity",color="black")+theme_classic() +
  scale_fill_manual(values = c(color.R, color.NR), breaks = c("nCRT-R", "nCRT-NR")) +
  geom_hline(yintercept= median(sampleTableCount.Clinical$MATH_score), linetype="dashed", color = "black")+
  annotate(geom="text", x=30, y=50, label=paste("median MATH_score= ",sprintf("%.2f",median(sampleTableCount.Clinical$MATH_score))), color="black")+
  xlab("SAMPLES") +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
gridExtra::grid.arrange(p, bottom=paste0("nCRT-R= ",aux[2], "; nCRT-NR= ",aux[1]))

pdf("output-Mutect2_SomaticFilters.2024-02-05/plots/barplot.math_score.pdf", width = 10, height = 6)
gridExtra::grid.arrange(p, bottom=paste0("MATH-score distribution from 52 samples\n nCRT-NR: ",aux[1],"; nCRT-R: ",aux[2] ))
dev.off()

aux<-(unique(somatic.52[,c("index", "ExonicFunc.refGene")]))  
aux<-as.data.frame(table(aux$ExonicFunc.refGene))
aux$perc <- (aux$Freq/sum(aux$Freq)) *100
aux

# TESTE DE NORMALIDADE DOS DADOS
shapiro.test(sampleTableCount.Clinical$MATH_score)
#==============================================================================#

# PLOT BOXPLOT MATH_score CORRELATION CLINICAL
#==============================================================================#
table(sampleTableCount.Clinical$ClinicalOutcome)
aux<-(table(sampleTableCount.Clinical$ClinicalOutcome))

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
pdf("output-Mutect2_SomaticFilters.2024-02-05/plots/boxplot.math_score.pdf", width = 4, height = 4)
gridExtra::grid.arrange(p, bottom= paste(result_test$method,"; p-value: ",sprintf("%.4f",result_test$p.value)))
dev.off()

## BARPLOT TMB ####
#==============================================================================#
p<- ggplot(data=sampleTableCount.Clinical, aes(x= reorder(ID_Exoma, -TMB), y=TMB, fill=ClinicalOutcome)) +
  geom_bar(stat="identity", color="black")+theme_classic() +
  scale_y_sqrt()+
  scale_fill_manual(values = c(color.R, color.NR), breaks = c("nCRT-R", "nCRT-NR")) +
  geom_hline(yintercept= median(sampleTableCount.Clinical$TMB), linetype="dashed", color = "black")+
  annotate(geom="text", x=30, y=50, label=paste("median TMB= ",sprintf("%.2f",median(sampleTableCount.Clinical$TMB))), color="black")+
  xlab("SAMPLES")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gridExtra::grid.arrange(p, bottom="TMB distribution from 52 samples")

pdf("output-Mutect2_SomaticFilters.2024-02-05/plots/barplot.tmb.pdf", width = 10, height = 4)
gridExtra::grid.arrange(p, bottom="TMB distribution from 52 samples")
dev.off()

sampleTableCount.Clinical$Perc_TMB <- sampleTableCount.Clinical$TMB/max(sampleTableCount.Clinical$TMB)

# TESTE DE NORMALIDADE DOS DADOS
shapiro.test(sampleTableCount.Clinical$TMB)
#==============================================================================#

# PLOT BOXPLOT TMB CORRELATION CLINICAL
#==============================================================================#
aux<-(table(sampleTableCount.Clinical$ClinicalOutcome))
p <- ggplot(sampleTableCount.Clinical  , aes(x=ClinicalOutcome, y=TMB, fill=ClinicalOutcome)) + 
  geom_boxplot()+ 
  scale_y_sqrt()+
  scale_fill_manual(values = c(color.R, color.NR), breaks = c("nCRT-R", "nCRT-NR")) +
  annotate("text", label = paste('N = ',aux[2], sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = ',aux[1], sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
   theme_classic()+theme(legend.position="none")
result_test <- wilcox.test(TMB ~ ClinicalOutcome, data = sampleTableCount.Clinical)
gridExtra::grid.arrange(p, bottom= paste(result_test$method,";\n p-value: ",sprintf("%.4f",result_test$p.value)))

pdf("output-Mutect2_SomaticFilters.2024-02-05/plots/boxplot.tmb.pdf", width = 4, height = 4)
gridExtra::grid.arrange(p, bottom= paste(result_test$method,";\n p-value: ",sprintf("%.4f",result_test$p.value)))
dev.off()

# sem ROP-107 e ROP-83
#==============================================================================#
aux<-(table(sampleTableCount.Clinical$ClinicalOutcome))
p <- ggplot(sampleTableCount.Clinical2 , aes(x=ClinicalOutcome, y=TMB, fill=ClinicalOutcome)) + 
  geom_boxplot()+ 
  # scale_y_sqrt()+
  scale_fill_manual(values = c(color.R, color.NR), breaks = c("nCRT-R", "nCRT-NR")) +
  annotate("text", label = paste('N = ',aux[2], sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = ',aux[1], sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
   theme_classic()+theme(legend.position="none")
result_test <- wilcox.test(TMB ~ ClinicalOutcome, data = sampleTableCount.Clinical2)
gridExtra::grid.arrange(p, bottom= paste(result_test$method,";\n p-value: ",sprintf("%.4f",result_test$p.value)))

pdf("output-Mutect2_SomaticFilters.2024-02-05/plots/boxplot.tmb.semROP-107-89.pdf", width = 4, height = 4)
gridExtra::grid.arrange(p, bottom= paste(result_test$method,";\n p-value: ",sprintf("%.4f",result_test$p.value)))
dev.off()

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

sp <- ggscatter(sampleTableCount.Clinical, x = "MATH_score", y = "TMB",
                add = "reg.line",  # Add regressin line
                cor.method = "spearman",
                color = "white",
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 30)+ 
  geom_point(aes(color = ClinicalOutcome))+
  theme(legend.position="top")
gridExtra::grid.arrange(sp, bottom="Spearman Test (52 samples)")

#sem outliers
sp <- ggscatter(sampleTableCount.Clinical2, x = "MATH_score", y = "TMB",
                add = "reg.line",  # Add regressin line
                cor.method = "spearman",
                color = "white",
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 40)+ 
  geom_point(aes(color = ClinicalOutcome))+
  theme(legend.position="top")
gridExtra::grid.arrange(sp, bottom="Spearman Test (50 samples) excluded outliers")

# RESPONSE 1
sp <- ggscatter(sampleTableCount.Clinical2, x = "MATH_score", y = "TMB",
                add = "reg.line",  # Add regressin line
                cor.method = "spearman",
                color = "white",
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman",  aes(color = ClinicalOutcome), label.x = 40)+
  geom_point(aes(color = ClinicalOutcome))+
  theme(legend.position="top")
gridExtra::grid.arrange(sp, bottom="Spearman Test (50 samples)")


#==============================================================================#
# TOP 50 GENES MAIS MUTADOS ####
#==============================================================================#
dgenes<-fread("../resultados_Mutect2-2023/input_Somatic-Filter/driver-genes.txt", header = F, col.names = "gene")
dMMR<-fread("../resultados_Mutect2-2023/input_Somatic-Filter/dMMR_genes.txt", header = F, col.names = "gene")

mut.genes=unique(dplyr::select(somatic.52, c("index","Gene.refGene")))
aux2 <- as.data.frame(table(mut.genes$Gene.refGene))
names(aux2) <- c("Gene", "Freq_Mutation")
aux2 <- (aux2[order(aux2$Freq, decreasing = T),])
aux2$class <-'none'
aux2$class[aux2$Gene %in% dgenes$gene] <-'driver'
aux2$class[aux2$Gene %in% dMMR$gene] <-'dMMR'
##fwrite(aux2, "output_Somatic-Filter/countMut_by_Genes.tsv", quote = F, sep="\t")
aux2=head(aux2,50)
top50<- aux2$Gene
p<-ggplot(data=aux2, aes(x= reorder(Gene, -Freq_Mutation), y=Freq_Mutation, fill=class)) +
  geom_bar(stat="identity" )+
  # scale_y_sqrt()+
  theme_classic() +xlab("Gene")
p +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#==============================================================================#
# TOP 50 GENES MAIS MUTADOS 50 samples ####
#==============================================================================#
dgenes<-fread("../resultados_Mutect2-2023/input_Somatic-Filter/driver-genes.txt", header = F, col.names = "gene")
dMMR<-fread("../resultados_Mutect2-2023/input_Somatic-Filter/dMMR_genes.txt", header = F, col.names = "gene")
somatic.50 <- subset(somatic.52, !(ID_Exoma %in% c("ROP-83", "ROP-107")), )

mut.genes=unique(dplyr::select(somatic.50, c("index","Gene.refGene")))
aux2 <- as.data.frame(table(mut.genes$Gene.refGene))
names(aux2) <- c("Gene", "Freq_Mutation")
aux2 <- (aux2[order(aux2$Freq, decreasing = T),])
aux2$class <-'none'
aux2$class[aux2$Gene %in% dgenes$gene] <-'driver'
aux2$class[aux2$Gene %in% dMMR$gene] <-'dMMR'
##fwrite(aux2, "output_Somatic-Filter/countMut_by_Genes.tsv", quote = F, sep="\t")
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
# #fwrite(target_file, "input_Somatic-Filter/xgen-exome-research-panel-v2-targets-hg38.cds", quote = F, col.names = T, row.names = F)
# #fwrite(as.data.frame(unique(target_file$gene.name)), "input_Somatic-Filter/xgen-exome-research-panel-v2-targets-hg38.list", quote = F, col.names = F, row.names = F)
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
# 52 samples
#==============================================================================#
all.mut.flt<-fread(paste0(output,"/mutation_codRegion_Freq005.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)
all.mut.flt.52<- subset(all.mut.flt, SAMPLE %in% unique(sampleTableCount.Clinical$SAMPLE))

# aux<-(mutations2[which(mutations2$Ref %in% c("A", "T","C", "G")),])
# aux<-(aux[which(aux$Alt %in% c("A", "T","C", "G")),])
#  unique(aux[,c('Ref', 'Alt')])

# targert_file<- fread("input_Somatic-Filter/xgen-exome-research-panel-v2-targets-hg38.bed ")
# targert_file <- separate(data = targert_file, col = V4, into = c("Gene","id"), sep = "_")
mutations <- all.mut.flt.52[,c("SAMPLE", "Chr", "Start", "Ref", "Alt")]
dndsout = dndscv(mutations, refdb="../resultados_Mutect2/input_Somatic-Filter/RefCDS_human_GRCh38.p12_dNdScv.0.1.0.rda", cv=NULL,
                 # gene_list = unique(targert_file$Gene),
                 max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf )
mutations2 <-mutations
mutations2$Chr <- gsub("chr", "", mutations2$Chr)
dndsout2 = dndscv(mutations2, refdb="../resultados_Mutect2/input_Somatic-Filter/RefCDS_human_GRCh38_GencodeV18_recommended.rda", cv=NULL,
                 # gene_list = unique(targert_file$Gene),
                 max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf )

mutations3 <- all.mut.flt.52[,c("SAMPLE", "CHR", "POS", "REF", "ALT")]
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

#fwrite(sel_cv[sel_cv$pglobal_cv<0.05, ], paste0(output,"/somatic.52.mut_signif_dndscv.tsv"), quote = F, sep="\t")
#fwrite(signif_genes, paste0(output,"/somatic.52.mut_signif_genes.tsv"), quote = F, sep="\t")

#==============================================================================#
#                ONCOPLOT  52 samples ####
#==============================================================================#
signif_genes <- fread(file=paste0(output,"/somatic.52.mut_signif_genes.tsv"))

somatic.maf<-fread(paste0(output,"/mutation_codRegion_Freq005.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)
somatic.maf.52<- subset(somatic.maf, SAMPLE %in% unique(sampleTableCount.Clinical$SAMPLE))
somatic.maf.52 <- merge(somatic.maf.52, sampleTableCount.Clinical[,c("SAMPLE","ID_Exoma")], by= "SAMPLE", all.x = T)
# somatic.maf<- subset(somatic.maf.52, select=c("CHR", "POS", "REF", "ALT", "Gene.refGene","ID_Exoma","ExonicFunc.refGene"))
somatic.maf<- subset(somatic.maf.52, select=c("Chr", "Start",  "End","Ref", "Alt", "Gene.refGene","ID_Exoma","ExonicFunc.refGene"))

somatic.maf<- somatic.maf %>% 
  mutate(Gene.refGene = strsplit(as.character(Gene.refGene), ";")) %>% 
  unnest(Gene.refGene)

# DESCOMENTE ABAIXO PARA PLOTAR SEM MUTAÇOES SINONIMAS
# somatic.maf<- subset(somatic.52, select=c("Chr", "Start",  "End","Ref", "Alt", "Gene.refGene","ID_Exoma","ExonicFunc.refGene"))
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

aux<-unique(subset(somatic.52, class == "dMMR",select = c("SAMPLE","class") ))
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
         clinicalFeatures = c("ClinicalOutcome"), 
         showTumorSampleBarcodes = T, SampleNamefontSize = 0.8,
         removeNonMutated = F,  drawColBar = T,
         fontSize = 0.6, annoBorderCol = "white",
         annotationFontSize = 1, 
         legend_height =2, 
         legendFontSize = 1)


pdf("output-Mutect2_SomaticFilters.2024-02-05/plots/oncoplot.pdf", width = 8, height = 5)
oncoplot(maf = maf, 
         genes = subset(signif_genes, qglobal_cv <= 0.1, select=gene_name)$gene_name,
         sortByAnnotation = T,
         clinicalFeatures = c("ClinicalOutcome"), 
         showTumorSampleBarcodes = T, SampleNamefontSize = 0.8,
         removeNonMutated = F,  drawColBar = T,
         fontSize = 0.6, annoBorderCol = "white",
         annotationFontSize = 1, 
         legend_height =2, 
         legendFontSize = 1)
dev.off()
#==============================================================================#
# DNDS_CV 2 - sem ROP-107 e ROP-83####
#==============================================================================#
all.mut.flt<-fread(paste0(output,"/mutation_codRegion_Freq005.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)
all.mut.flt.52<- subset(all.mut.flt, SAMPLE %in% unique(sampleTableCount.Clinical$SAMPLE))
all.mut.flt.50<-subset(all.mut.flt.52, 
                        !(SAMPLE %in% c("ROP-83-ExC85-xgenV2_S54", "ROP-107-ExC85-xgenV2_S71")),)
mutations2 <- all.mut.flt.50[,c("SAMPLE", "Chr", "Start", "Ref", "Alt")]
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
#fwrite(sel_cv[sel_cv$pglobal_cv<0.05, ], paste0(output,"/somatic.50.mut_signif_dndscv.tsv"), quote = F, sep="\t")
#fwrite(signif_genes, paste0(output,"/somatic.50.mut_signif_genes.tsv"), quote = F, sep="\t")

#==============================================================================#
#                ONCOPLOT  50 samples ####
#==============================================================================#
signif_genes <- fread(file=paste0(output,"/somatic.50.mut_signif_genes.tsv"))

somatic.maf<-fread(paste0(output,"/mutation_codRegion_Freq005.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)
somatic.maf.50<- subset(somatic.maf, SAMPLE %in% unique(sampleTableCount.Clinical2$SAMPLE))
somatic.maf.50 <- merge(somatic.maf, sampleTableCount.Clinical2[,c("SAMPLE","ID_Exoma")], by= "SAMPLE", all.x = T)
somatic.maf<- subset(somatic.maf.50, select=c("Chr", "Start",  "End","Ref", "Alt", "Gene.refGene","ID_Exoma","ExonicFunc.refGene"))

somatic.maf<- somatic.maf %>% 
  mutate(Gene.refGene = strsplit(as.character(Gene.refGene), ";")) %>% 
  unnest(Gene.refGene)

# DESCOMENTE ABAIXO PARA PLOTAR SEM MUTAÇOES SINONIMAS
# somatic.maf<- subset(somatic.50, select=c("Chr", "Start",  "End","Ref", "Alt", "Gene.refGene","ID_Exoma","ExonicFunc.refGene"))
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

aux<-unique(subset(somatic.50, class == "dMMR",select = c("SAMPLE","class") ))
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
         clinicalFeatures = c("ClinicalOutcome"), 
         showTumorSampleBarcodes = T, SampleNamefontSize = 0.8,
         removeNonMutated = F,  drawColBar = T,
         fontSize = 0.6, annoBorderCol = "white",
         annotationFontSize = 1, 
         legend_height =2, 
         legendFontSize = 1)


#==============================================================================#
#  pheatmap - GENES with dMMR  
#==============================================================================#
mut.genes=unique(dplyr::select(somatic.52, c("Gene.refGene","ID_Exoma", "class")))
aux<- as.data.frame(unique(subset(mut.genes, select= ID_Exoma)))

mut.genes.heatmap= as.data.frame(table(unique(dplyr::select(somatic.52, c("Gene.refGene","ID_Exoma")))))
mut.genes.heatmap<- mut.genes.heatmap[mut.genes.heatmap$Gene.refGene %in% dMMR$gene,]
mut.genes.heatmap<- mut.genes.heatmap[mut.genes.heatmap$ID_Exoma %in% aux$ID_Exoma,]
df<- as.matrix(reshape2::dcast(mut.genes.heatmap, Gene.refGene ~ ID_Exoma,  value.var = "Freq"))
rownames(df) <- df[,1]
df<-type.convert(df[,-1], as.is = T)

## Juntando tabelas TMB_MATH_#MUT +  dados_clinicos 
my_sample_col <- sampleTableCount.Clinical
my_sample_col <- data.frame(subset(my_sample_col, select = c("ID_Exoma","ClinicalOutcome")))

row.names(my_sample_col) <- my_sample_col[,1]
my_sample_col<-my_sample_col[order(my_sample_col$ClinicalOutcome ),]
df_ordered <- df[,rownames(my_sample_col)]
my_sample_col<-subset(my_sample_col, select=-1)

pheatmap((df_ordered), color=colorRampPalette(c("white", "darkred"))(100),
         annotation_col = my_sample_col,
         border_color = "grey60",
         cluster_rows= F, cluster_cols= F, fontsize = 8,
         legend = F )


#==============================================================================#
#  pheatmap - GENES with dMMR  67samples
#==============================================================================#
# mut.genes=unique(select(somatic.50, c("Gene.refGene","ID_Exoma", "class")))
# aux<- as.data.frame(unique(subset(mut.genes, select= ID_Exoma)))
# 
# mut.genes.heatmap= as.data.frame(table(unique(select(somatic.50, c("Gene.refGene","ID_Exoma")))))
# mut.genes.heatmap<- mut.genes.heatmap[mut.genes.heatmap$Gene.refGene %in% dMMR$gene,]
# mut.genes.heatmap<- mut.genes.heatmap[mut.genes.heatmap$ID_Exoma %in% aux$ID_Exoma,]
# df<- as.matrix(reshape2::dcast(mut.genes.heatmap, Gene.refGene ~ ID_Exoma,  value.var = "Freq"))
# rownames(df) <- df[,1]
# df<-type.convert(df[,-1], as.is = T)
# 
# my_sample_col <- sampleTableCount.Clinical2
# 
# my_sample_col <- data.frame(subset(my_sample_col, select = c("ID_Exoma","ClinicalOutcome","Response2","Response3")))
# row.names(my_sample_col) <- my_sample_col[,1]
# my_sample_col<-subset(my_sample_col, select=-1)
# 
# my_sample_col<-my_sample_col[order(my_sample_col$ClinicalOutcome ),]
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
# SALVAR ARQUIVOS BED PARA FILTRAR VCF ####
#==============================================================================#
samples<-unique(sampleTableCount.Clinical$SAMPLE)
dir.create(paste0(output,"mutations_to_filter_fromVCF/"))

for (sp in samples) {
  BED<-subset(somatic.52, SAMPLE==sp, select = c("CHR", "POS","REF","ALT","index"))
  BED<-unique(BED)
  #fwrite(BED,paste0(output,"mutations_to_filter_fromVCF/",sp,".bed"), sep = "\t", col.names = T, row.names = F, quote = F)
}



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
mutados<- unique(subset(somatic.52, class == "dMMR", select = SAMPLE))
final_Sig[final_Sig$SAMPLE %in% mutados$SAMPLE,"dMMR"] <-"dMMR"
# final_Sig<-merge(final_Sig, subset(somatic_Samples.resp, select = c("SAMPLE","ClinicalOutcome","Response3")), by="sample", all.y = T)
final_Sig<-merge(final_Sig, subset(sampleTableCount.Clinical, select = c("SAMPLE","ClinicalOutcome","Response3")), by="SAMPLE", all.y = T)
final_Sig$perc <-as.numeric(final_Sig$perc)
final_Sig$cosmic_sig <- "sig_others"
final_Sig[final_Sig$sig_id %in% c("6","15","20","26"),"cosmic_sig"] <-"MMR_sig"

#fwrite(final_Sig, "table_signatures.tsv", sep = "\t", quote = F, col.names = T, row.names = F)

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

aux<-unique(subset(somatic.52, class == "dMMR",select = c("SAMPLE","class") ))
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

score_MMR_sig <- subset(mutsig, cosmic_sig=="MMR_sig") %>% group_by(Sample, ClinicalOutcome, Response3) %>% 
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
#ClinicalOutcome vs MMR_sig
p<-ggplot(aux[aux$score_MMR_sig >0, ], 
          aes(x=ClinicalOutcome, y=score_MMR_sig)) + 
  geom_boxplot(outlier.shape = NA) + stat_compare_means() +
  annotate("text", label = paste('N = 10', sep = ''), x = 2, y = 0, angle = 0, size = 3.75) +
  annotate("text", label = paste('N = 21', sep = ''), x = 1, y = 0, angle = 0, size = 3.75) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_classic() +theme(legend.position="left")
gridExtra::grid.arrange(p, bottom="")
print(table(aux[aux$score_MMR_sig >0, ]$ClinicalOutcome))
table(aux$ClinicalOutcome)

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
#somatic.52<-as.data.frame(fread(paste0(output,"/somatic_mutation.72samples.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))
# 
# aux1 <- dplyr::select((somatic.52  ), SAMPLE, MAF)
# aux2 <- as.data.frame(table(aux1), stringsAsFactors = F)
# aux2$MAF <- as.double(aux2$MAF)
# aux2 <- merge(aux2, sampleTableCount.Clinical[,c("SAMPLE","ID_Exoma","ClinicalOutcome")], by= "SAMPLE", all.x = T)

# MAFs das amostras outliers    ####
#==============================================================================#
# # nCRT-R ClinicalOutcome
# P<- ggplot(data=aux2,
#            aes(x=MAF, y=Freq, group = ID_Exoma, colour= ClinicalOutcome))+
#   geom_col(colour= color.R) +
#   # scale_y_sqrt()+
#   scale_x_continuous(breaks=seq(0,1,0.1))+
#   facet_wrap(vars(ID_Exoma))+theme_minimal()+
#   theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")
# P

#==============================================================================#
# Fisher Test  - Drivers genes, dMMR ####
#==============================================================================#
classeGenes <- c("DriverGene", "dMMR")
respostas <- c("ClinicalOutcome")

for (geneClass in classeGenes) {
  for (resp in respostas) {
    # resp <- "ClinicalOutcome"
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

dndsGenes <-fread(paste0(output,"somatic.50.mut_signif_genes.tsv"), header = T, na.strings=c("NA"), fill=TRUE, check.names = FALSE)
# dndsGenes <- subset(dndsGenes, qglobal_cv <= 0.1, select = "gene_name")

for (resp in respostas) {
  for (gene in dndsGenes$gene_name) {
    # resp <- "ClinicalOutcome"
    # gene="APC"
    aux <-sampleTableCount.Clinical
    aux$mutGene <-"not_dNdSGene"
    mutados<- unique(subset(somatic.52,  Gene.refGene == gene, select = ID_Exoma))
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
respostas <- c("ClinicalOutcome", "Response2", "Response3")

for (id in List_Ass) {
  for (resp in respostas) {
  # id="28"
  # resp= "ClinicalOutcome"
  aux<-as.data.frame(unique(subset(final_Ass,  sig_id==id, select=c("SAMPLE"))), stringsAsFactors = F)
  table_sig$MMR_sig = "0"
  table_sig[table_sig$SAMPLE %in% aux$SAMPLE]$MMR_sig <-"1"

  # dt3<- matrix(table(table_sig[,c("MMR_sig", "ClinicalOutcome")]), nr=2)
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
respostas <- c("ClinicalOutcome", "Response2", "Response3")

final_Ass<-subset(final_Ass, dMMR == "dMMR")
for (id in List_Ass) {
  for (resp in respostas) {
    # id="4"
    # resp= "Response2"
    aux<-as.data.frame(unique(subset(final_Ass,  sig_id==id, select=c("SAMPLE"))), stringsAsFactors = F)
    table_sig$MMR_sig = "0"
    table_sig[table_sig$SAMPLE %in% aux$SAMPLE]$MMR_sig <-"1"
    
    # dt3<- matrix(table(table_sig[,c("MMR_sig", "ClinicalOutcome")]), nr=2)
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
# INSERIR INFO DRIVER GENES MUTADOS E pMMR  ####
#==============================================================================#

sampleTableCount.Clinical$DriverGene <-"notDriver"
mutados<- unique(subset(somatic.52, class == "driver", select = ID_Exoma))
sampleTableCount.Clinical[sampleTableCount.Clinical$ID_Exoma %in% mutados$ID_Exoma]$DriverGene <-"RectalDriver"

sampleTableCount.Clinical$dMMR <-"pMMR"
mutados<- unique(subset(somatic.52, class == "dMMR", select = ID_Exoma))
sampleTableCount.Clinical[sampleTableCount.Clinical$ID_Exoma %in% mutados$ID_Exoma]$dMMR <-"dMMR"

#==============================================================================#
# CARREGANDO INFORMAÇÕES DE SVs ####
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
respostas <- c("ClinicalOutcome")

for (sv in eventos_SV) {
  for (resp in respostas) {
    # sv="status_SV"
    # resp= "ClinicalOutcome"
    
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


#==============================================================================#
# ANALISE MULTIVARIADA ####
#==============================================================================#

immuno<-fread(paste0(output,"/ROP_immuno_69pat_final.tsv"), quote = F, sep="\t")
immuno<-subset(immuno, Sample_Name %in% sampleTableCount.Clinical$ID_Exoma, 
               select = c('Sample_Name','Clones','Heterozygosis','HED_Score','Mean_class','NAL','nal_bin','NAL_snvs','NAL_indel','HLA.mut'))

dados<-sampleTableCount.Clinical
dados$ClinicalOutcome <-ifelse(dados$ClinicalOutcome == "nCRT-NR", 1, 0) # USAR VITAL.STATUS OU RESPOSTA, RESPOSTA2
dados$Distance.from.anal.verge <- gsub("cm", "", dados$Distance.from.anal.verge)
dados$Tumor.size <- gsub("cm", "", dados$Tumor.size)
dados$Tumor.size <- gsub("xx", NA, dados$Tumor.size)
dados$Distance.from.anal.verge <- gsub("xx", NA, dados$Distance.from.anal.verge)
dados$CEA.inicial <- gsub("xx", NA, dados$CEA.inicial)
dados$Tratamento <- gsub("normal", "nCRT", dados$Tratamento)
dados$Tratamento <- gsub("extended", "TNT", dados$Tratamento)
dados$Distance.from.anal.verge <-as.numeric(dados$Distance.from.anal.verge)
dados<- merge(dados, immuno, by.x = "ID_Exoma", by.y = "Sample_Name")
#MATH-score pROC
dados$MATH_score.= "Low_MATH"
colnames(dados)[colnames(dados) == "Tratamento"] <- "Treatment.Protocol"

dados$MATH_score.= ifelse(dados$MATH_score >= 30.75718, "High_MATH", dados$MATH_score.)
dados$Clones <- gsub("-", NA, dados$Clones)

#==============================================================================#
# ANALISE MULTIVARIADA #### 
# Continua no script "Analise multivariada.R"
#==============================================================================#


#==============================================================================#
# EXEMPLO MAIS SIMPLE DE CURVA ROC ####
#==============================================================================#
aux<-as.data.frame(sampleTableCount.Clinical[, c('MATH_score','ClinicalOutcome')],)
aux$ClinicalOutcome <-ifelse(aux$ClinicalOutcome == "nCRT-NR", 1, 0) # USAR VITAL.STATUS OU RESPOSTA, RESPOSTA2
dados<-aux

# Ajuste a curva ROC
roc_curve <- roc(dados$ClinicalOutcome, dados$MATH_score)
# Plote a curva ROC
plot(roc_curve, main = "Curva ROC", col = "blue")

# Adicione informações ao gráfico (opcional)
text(0.5, 0.5, paste("AUC =", round(auc(roc_curve), 2)), adj = c(0.5, -0.5), col = "black", cex = 1.2)

# Get the best threshold
coords(roc_curve, "best", ret="threshold", transpose = FALSE)
coords(roc_curve, "best", ret=c("threshold", "specificity", "sensitivity", "accuracy",
                                "precision", "recall"), transpose = T)


# =============================================================================#
#                        CURVA ROC  TMB e MATH_score - RESPOSTA ####
# =============================================================================#
library(pROC)
curvaROC <- function(tabela, coluna){
  colunas<-c(coluna,"ClinicalOutcome")
  # Converte para data.table
  setDT(tabela)
  
  # Verifica se as colunas existem no data.table
  colunas_existentes <- colunas[colunas %in% colnames(tabela)]
  
  if (length(colunas_existentes) > 0) {
    # Cria um novo data.table com as colunas desejadas
    novo_data_table <- tabela[, ..colunas_existentes, with = FALSE]
   
    aux<-as.data.frame(novo_data_table)

    a<-ggplot(aux, aes_string(x = coluna, fill = "ClinicalOutcome")) +
      #geom_density(position = "identity", alpha = 0.5) +
      geom_histogram(position = "identity", alpha = 0.5, bins = 20) +
      scale_fill_manual(values = c("red", "blue")) +
      ggtitle("Distribuition (52 samples)") + theme_classic()+
      xlab(coluna) #+ ylab("Freq.")
    
    # ============================#
    # curva roc Completa vs Incompleta
    # ============================#
    aux$ClinicalOutcome <-ifelse(aux$ClinicalOutcome == "nCRT-NR", 1, 0) # USAR VITAL.STATUS OU RESPOSTA, RESPOSTA2
    
    #=============================#
    # Dividindo o dataset em treino e teste
    # Criar o subset de treino
    train <- aux %>% dplyr::sample_frac(.70)
    # Criar o subset de teste com antijoin (pega tudo que não pertence)
    test <- aux 
    
    # Rodando o modelo
    formula <- as.formula(paste("ClinicalOutcome ~", coluna))
    fit_reg_log <- 
      glm(formula,
          family = binomial(link = 'logit'),
          data = train)
    
    # Fazendo as predições
    pred_reg_log <- predict(fit_reg_log, newdata = test, type = "response")
    
    # Organizando a tabela de dados para calcular as métricas da curva ROC
    pred_roc_reg_log <- 
      dplyr::tibble(
        pred_reg_log,
        "survived" = as.factor(as.numeric(test$ClinicalOutcome)-1)
      ) %>% arrange(desc(pred_reg_log))
    
    # Criando objeto com as métricas para curva ROC
    roc_reg_log <- pROC::roc(pred_roc_reg_log$survived , pred_roc_reg_log$pred_reg_log, percent = TRUE)
    # Se desejar, é possível (e bem simples) utilizar o próprio pacote pROC para plotar a curva ROC.
    par(pty = "s")
    plot.roc(
      main = "ROC Curve (52 samples) \n nCRT-R vs nCRT-NR",
      roc_reg_log,
      print.auc = TRUE,
      legacy.axes = TRUE
      # xlab = "Taxa de Falso Positivo (100 - Especificidade)",
      # ylab = "Taxa de Verdadeiro Positivo (Sensibilidade)"
    )
    
    # Retorna o novo data.table se você quiser usá-lo fora da função
    return(a)
  } else {
    cat("Nenhuma coluna especificada existe no data.table.\n")
    return(NULL)
  }
}
curvaROC(sampleTableCount.Clinical, c("MATH_score") )
#curvaROC(sampleTableCount.Clinical2, c("TMB") )



#==============================================================================#
# Fisher Test  - ClinicalOutcome vs Treatment.Protocol ####
#==============================================================================#
colnames(dados)[colnames(dados) == "Tratamento"] <- "Treatment.Protocol"

col.vars<- "MATH_score."
for (col.var in col.vars) {
  resp <- "ClinicalOutcome"
  # col.var<- "dMMR"
  dt3<- matrix(table(subset(dados, select = c(col.var, resp))), nr=2)
  if(ncol(dt3)>1){
    resultado<-fisher.test(dt3)
    valor_p <- resultado$p.value
    if(valor_p <= 5){
      print("=============================================================")
      print(col.var)
      print(resp)
      print(as.matrix(table(subset(dados, select = c(col.var, resp))), nr=2))
      print(resultado)
      print("=============================================================")
      print(chisq.test(dt3))
      print("=============================================================")
    }
  }
}

