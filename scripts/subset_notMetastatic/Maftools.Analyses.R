library(maftools)
library(data.table)
library(dplyr)
library(ggplot2)
library('pheatmap')
library('NMF')
library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)

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
maf.Tab <- read.maf(maf = syn.Mutation.MAF, 
                    clinicalData = Clinical)

pdf(paste0(output,"Figures/summary.MAFTOOLS.pdf"), width =14, height = 8)
plotmafSummary(maf = maf.Tab, rmOutlier = TRUE, addStat = 'median', showBarcodes = T,
               dashboard = TRUE, titvRaw = F , textSize = 0.6)
maf.Tab.titv = titv(maf = maf.Tab, plot = FALSE, useSyn = F)
#plot titv summary
plotTiTv(res = maf.Tab.titv,showBarcodes = T)


# Generate a bar plot displaying distinct ExonicFunc mutation types
aux<-as.data.frame(unique(syn.Mutation.MAF[,c("index", "ExonicFunc.refGene")]))
aux <- as.data.frame(table(aux$ExonicFunc.refGene))
aux <- aux %>%
  arrange(desc(Freq)) %>%
  mutate(Var1 = factor(Var1, levels = unique(Var1)))
colnames(aux) <- c("ExonicFunc", "count")
aux <- aux %>%  mutate(proportion = count / sum(count) * 100)

p<-ggplot(data = aux, aes(x = ExonicFunc, y = count, fill = ExonicFunc)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = paste(count, sprintf("(%.1f%%)", proportion)), y = 10000),
            hjust = 0.1, color = "black") +
  coord_flip() +  theme_minimal() +
  labs(title = "30.282 Somatic Variants",
       x = "ExonicFunc Mutation Type",
       y = "Count") +
  scale_fill_grey(start = 0.8, end = 0.2) +
  theme(legend.position = "none")
gridExtra::grid.arrange(p, bottom="ExonicFunc.refGene distribution")

# Generate a bar plot displaying distinct Variant_Classification
aux<-as.data.frame(unique(syn.Mutation.MAF[,c("index", "Variant_Classification")]))
aux <- as.data.frame(table(aux$Variant_Classification))
aux <- aux %>%
  arrange(desc(Freq)) %>%
  mutate(Var1 = factor(Var1, levels = unique(Var1)))
colnames(aux) <- c("Variant_Classification", "count")
aux <- aux %>%  mutate(proportion = count / sum(count) * 100)

p<-ggplot(data = aux, aes(x = Variant_Classification, y = count, fill = Variant_Classification)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = paste(count, sprintf("(%.1f%%)", proportion)), y = 10000),
            hjust = 0.1, color = "black") +
  coord_flip() +  theme_minimal() +
  labs(title = "30.282 Somatic Variants",
       x = "Variant_Classification",
       y = "Count") +
  scale_fill_grey(start = 0.8, end = 0.2) +
  theme(legend.position = "none")
gridExtra::grid.arrange(p, bottom="ExonicFunc.refGene distribution")



Mb <- 34126765/1000000 
TMB.maftools <-tmb(maf.Tab, captureSize = Mb, logScale = TRUE)
plot.new()
plot.new()
rop.mutload<- tcgaCompare(maf = maf.Tab, cohortName = 'ROP-SAMPLES',
                          logscale = TRUE, capture_size = Mb)
plot.new()
# 
# #grid.newpage()
# plot(rop.mutload = tcgaCompare(maf = subsetMaf(maf.Tab, tsb = samplesNotHiperMut$EXOME_ID),
#                           cohortName = 'ROP-SAMPLES', rm_hyper =T,
#                           logscale = TRUE, capture_size = Mb))

##9.2 Detecting cancer driver genes based on positional clustering
# maftools has a function oncodrive which identifies cancer genes (driver) from a given MAF. 
# oncodrive is a based on algorithm oncodriveCLUST which was originally implemented in Python. 
# Concept is based on the fact that most of the variants in cancer causing genes are enriched at few specific loci (aka hot-spots). 
# This method takes advantage of such positions to identify cancer genes. 
# If you use this function, please cite OncodriveCLUST article 7.
rop.sig = oncodrive(maf = maf.Tab, AACol = 'aaChange', minMut = 5, pvalMethod = 'zscore',
                    ignoreGenes = c("TTN"))
fwrite(rop.sig, paste0(output,"/rop.sig.tsv"), quote = F, sep="\t")

plotOncodrive(res = rop.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)
rop.sig[rop.sig$fdr <= 0.1, ]$Hugo_Symbol
lollipopPlot(maf = maf.Tab, gene = "ABCF1", AACol = "aaChange" )
lollipopPlot(maf = maf.Tab, gene = "OR4A15", AACol = "aaChange" )
lollipopPlot(maf = maf.Tab, gene = "CDH5", AACol = "aaChange" )
lollipopPlot(maf = maf.Tab, gene = "KRAS", AACol = "aaChange" )

## 9.3 Adding and summarizing pfam domains
# maftools comes with the function pfamDomains, which adds pfam domain information to the amino acid changes. 
# pfamDomain also summarizes amino acid changes according to the domains that are affected.
# This serves the purpose of knowing what domain in given cancer cohort, is most frequently affected. 
# This function is inspired from Pfam annotation module from MuSic tool 8.
rop.pfam = pfamDomains(maf = maf.Tab, AACol = 'aaChange', top = 10, baseName ="output-SomaticFilter/rop-pfam" )
#Protein summary (Printing first 7 columns for display convenience)
rop.pfam$proteinSummary[,1:7, with = FALSE]

dev.off()


#==============================================================================#
# DISTRIBUIÇÃO DAS VAFs   ####
#==============================================================================#

aux1 <- select(subset(syn.Mutation.MAF), SAMPLE, VAF)
aux2 <- as.data.frame(table(aux1), stringsAsFactors = F)
aux2$VAF <- as.double(aux2$VAF)
aux2 <- merge(aux2, Clinical[,c("Sample","EXOME_ID","Response1","Response2")], by.x= "SAMPLE", by.y = "Sample", all.x = T)
aux2 <- aux2[aux2$Freq >0,]

# Função para criar gráficos
create_plots <- function(data, response, color, output) {
  subset_data <- subset(data, Response1 == response)
  # Gráfico de barras
  p1 <- ggplot(subset_data, aes(x = VAF, y = Freq, group = EXOME_ID)) +
    geom_col(colour = color) +
    scale_y_sqrt() +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    facet_wrap(vars(EXOME_ID)) +
    theme_minimal() +
    labs(x = "VAF", y = "Density", title = paste("Density Plot of VAF -", response)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "top")
  
  # Gráfico de densidade
  p2 <- ggplot(subset_data, aes(x = VAF, group = EXOME_ID)) +
    geom_density(fill = color, alpha = 0.5) +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    facet_wrap(vars(EXOME_ID)) +
    theme_minimal() +
    labs(x = "VAF", y = "Density", title = paste("Density Plot of VAF -", response)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "top")
  
  # Boxplot
  p3 <- ggplot(subset_data, aes(x = VAF, y = "", fill = Response1)) +
    geom_boxplot(alpha = 0.5, fill = "gray") +
    geom_vline(xintercept = 0.2, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    facet_wrap(vars(EXOME_ID)) +
    labs(title = paste("Boxplot of VAF for each EXOME_ID -", response), x = "VAF", y = "VAF") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "top")
  
  # Salvar os gráficos em um arquivo PDF
  pdf(file = paste0(output, "Figures/VAF_Distribution_", gsub(" ", "_", response), ".synonimous.pdf"), width = 12, height = 6)
  print(p1)
  print(p2)
  print(p3)
  dev.off()
}

# Gerar gráficos para cada grupo de Response2
create_plots(aux2, "nCRT-R", "#08415c", output)
create_plots(aux2, "nCRT-NR", "#ef767a", output)


#==============================================================================#
# DNDS_CV ####
#==============================================================================#
library(dndscv)

#==============================================================================#
aux <- syn.Mutation.MAF[,c("SAMPLE", "CHR", "Start_Position","REF", "ALT")]
aux$CHR <- gsub("chr", "", aux$CHR)
dndsout = dndscv(aux, refdb="D:/PROJETOS-HSL_BP/resultados_Mutect2-2023/input_Somatic-Filter/RefCDS_human_GRCh38_GencodeV18_recommended.rda", 
                 cv=NULL, 
                 # gene_list = unique(targert_file$Gene),
                 max_coding_muts_per_sample = Inf, 
                 max_muts_per_gene_per_sample = Inf )
print(dndsout$globaldnds)
# In particular, very low estimates of θ (the overdispersion parameter), 
# particularly θ<1, may reflect problems with the suitability of the dNdScv model for the dataset.
print(dndsout$nbreg$theta)

signif_genes_localmodel = as.vector(dndsout$sel_loc$gene_name[dndsout$sel_loc$qall_loc<0.01])
# View(dndsout$sel_cv) # This is shown as an example but these results based on a few genes should not be trusted
sel_cv = dndsout$sel_cv
signif_genes = sel_cv[sel_cv$pglobal_cv<0.05, c("gene_name", "pglobal_cv" ,"qglobal_cv","qallsubs_cv")]
nrow(subset(signif_genes,pglobal_cv < 0.05))
nrow(subset(signif_genes,qglobal_cv < 0.1))
nrow(subset(signif_genes,qglobal_cv < 0.05))
nrow(subset(signif_genes,qallsubs_cv < 0.05)) # avaliando sem os Indels
fwrite(sel_cv[sel_cv$pglobal_cv< 0.05, ], paste0(output,"/dndscv_pvalue_0.05.tsv"), quote = F, sep="\t")
fwrite(sel_cv[sel_cv$qglobal_cv< 0.1, ], paste0(output,"/dndscv_FDR_0.1.tsv"), quote = F, sep="\t")


# MutSigCV  ####
#==============================================================================#
# Prepare MAF file for MutSigCV analysis
mutsig.corrected = prepareMutSig(maf = maf.Tab)


#==============================================================================#
# dNds_genes PLOTS ####
#==============================================================================#
dnds_genes <- fread(file=paste0(output,"/dndscv_FDR_0.1.tsv"))

# ONCOPLOT  ####
#==============================================================================#
pdf(file = paste0(output,"Figures/dnds_genes.mutations.pdf"), width = 12, height = 6)
oncoplot(maf = maf.Tab, 
         genes = subset(signif_genes, qglobal_cv <= 0.05, select=gene_name)$gene_name,
         sortByAnnotation = T, 
#         pathways = "smgbp",collapsePathway = TRUE, topPathways = 10,gene_mar = 8,
         clinicalFeatures = c("Response1"), 
         showTumorSampleBarcodes = T, SampleNamefontSize = 0.8,
         removeNonMutated = F,  drawColBar = T,
         fontSize = 0.6, annoBorderCol = "white",
         annotationFontSize = 1, 
         legend_height =2, 
         legendFontSize = 1)

# dnds_genes Mutualmente/Exclusivamente mutados  ####
#==============================================================================#
somaticInteractions(maf = maf.Tab, genes = dnds_genes$gene_name, 
                    leftMar = 5, topMar = 5, 
                    pvalue = c(0.05, 0.05),fontSize = 0.6 )

# Looliplot dnds genes  ####
#==============================================================================#
dnds_genes$gene_name
lollipopPlot(maf = maf.Tab, gene = "KRAS", AACol = "aaChange" )
lollipopPlot(maf = maf.Tab, gene = "TP53", AACol = "aaChange" )
lollipopPlot(maf = maf.Tab, gene = "APC", AACol = "aaChange" )
lollipopPlot(maf = maf.Tab, gene = "FBXW7", AACol = "aaChange" )
lollipopPlot(maf = maf.Tab, gene = "SMAD4", AACol = "aaChange" )
lollipopPlot(maf = maf.Tab, gene = "ARID2", AACol = "aaChange" )
lollipopPlot(maf = maf.Tab, gene = "ZFP36L2", AACol = "aaChange" )
lollipopPlot(maf = maf.Tab, gene = "CD58", AACol = "aaChange" )

#9.6 Clinical enrichment analysis
# Performs fishers test on 2x2 contingenc
#==============================================================================#
rop.ce = clinicalEnrichment(maf = maf.Tab,  minMut = 5, 
                             clinicalFeature = 'Response1')
head(rop.ce$groupwise_comparision[p_value < 0.05])
plotEnrichmentResults(enrich_res = rop.ce, pVal = 0.05,  
                      geneFontSize = 0.5, annoFontSize = 0.6)
dev.off()

#==============================================================================#
# dNds_genes PLOTS 2 ####
#==============================================================================#
maf.dnds <- subsetMaf(maf.Tab, genes = subset(signif_genes, qglobal_cv <= 0.05, select=gene_name)$gene_name)

# dnds_genes Mutualmente/Exclusivamente mutados  ####
#==============================================================================#
pdf(file = paste0(output,"Figures/dnds_genes_2.mutations.pdf"), width = 12, height = 6)
oncoplot(maf = maf.dnds,
         sortByAnnotation = T , 
         #         pathways = "smgbp",collapsePathway = TRUE, topPathways = 10,gene_mar = 8,
         clinicalFeatures = c("Response1"), 
         showTumorSampleBarcodes = T, SampleNamefontSize = 0.8,
         removeNonMutated = F,  drawColBar = T,
         fontSize = 0.6, annoBorderCol = "white",
         annotationFontSize = 1, 
         legend_height =2, 
         legendFontSize = 1)

# dnds_genes Mutualmente/Exclusivamente mutados  ####
#==============================================================================#
somaticInteractions(maf = maf.dnds,    leftMar = 5, topMar = 5, 
                    pvalue = c(0.05, 0.05), fontSize = 0.6 )

# Looliplot dnds genes  ####
#==============================================================================#
dnds_genes$gene_name
lollipopPlot(maf = maf.dnds, gene = "KRAS", AACol = "aaChange" )
lollipopPlot(maf = maf.dnds, gene = "TP53", AACol = "aaChange" )
lollipopPlot(maf = maf.dnds, gene = "APC", AACol = "aaChange" )
lollipopPlot(maf = maf.dnds, gene = "FBXW7", AACol = "aaChange" )
lollipopPlot(maf = maf.dnds, gene = "SMAD4", AACol = "aaChange" )
lollipopPlot(maf = maf.dnds, gene = "ARID2", AACol = "aaChange" )
lollipopPlot(maf = maf.dnds, gene = "ZFP36L2", AACol = "aaChange" )
lollipopPlot(maf = maf.dnds, gene = "CD58", AACol = "aaChange" )

#9.6 Clinical enrichment analysis
# Performs fishers test on 2x2 contingenc
#==============================================================================#
dnds.ce = clinicalEnrichment(maf = maf.dnds,  minMut = 1, 
                            clinicalFeature = 'Response1')
head(dnds.ce$groupwise_comparision[p_value < 0.05])
plotEnrichmentResults(enrich_res = dnds.ce, pVal = 0.05, 
                      geneFontSize = 0.5, annoFontSize = 0.6)
dev.off()

#==============================================================================#
#9.5 Comparing two cohorts (MAFs)
#==============================================================================#
samples.nCRT_R <- Clinical[Clinical$Response1 == "nCRT-R",]$EXOME_ID
samples.nCRT_NR <- Clinical[Clinical$Response1 == "nCRT-NR",]$EXOME_ID
maf.nCRT_R<-subsetMaf(maf.Tab, tsb = samples.nCRT_R)
maf.nCRT_NR<-subsetMaf(maf.Tab, tsb = samples.nCRT_NR)

#Considering only genes which are mutated in at-least in 5 samples in one of the cohort to avoid bias due to genes mutated in single sample.
pt.vs.rt <- mafCompare(m1 = maf.nCRT_NR, m2 = maf.nCRT_R, pseudoCount = T,
                       m1Name = 'nCRT Not-Responders', m2Name = 'nCRT Responders',
                       minMut = 5)

pdf(file = paste0(output,"Figures/Comparing.Response.pdf"), width = 12, height = 6)
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05, fdr = 0.1)
genes.pval = pt.vs.rt$results[pt.vs.rt$results$pval <= 0.05, "Hugo_Symbol"]
genes.fdr = pt.vs.rt$results[pt.vs.rt$results$adjPval <= 0.1, "Hugo_Symbol"]

coOncoplot(m1 = maf.nCRT_NR, m2 = maf.nCRT_R, 
           m1Name = 'nCRT Not-Responders', m2Name = 'nCRT Responders', 
           genes = genes.fdr$Hugo_Symbol, removeNonMutated = TRUE)
coBarplot(m1 = maf.nCRT_NR, m2 = maf.nCRT_R, 
          m1Name = 'nCRT Not-Responders', m2Name = 'nCRT Responders', 
          genes = genes.fdr$Hugo_Symbol )
dev.off()

#==============================================================================#
#9.6 Clinical enrichment analysis
# Performs fishers test on 2x2 contingenc
#==============================================================================#
pdf(file = paste0(output,"Figures/Enrichment_Genes.Response.pdf"), width = 12, height = 6)
rop.ce = clinicalEnrichment(maf = maf.Tab, clinicalFeature = 'Response1')
rop.ce$groupwise_comparision[p_value < 0.05]
plotEnrichmentResults(enrich_res = rop.ce, pVal = 0.05, 
                      geneFontSize = 0.5, annoFontSize = 0.6)

# Mutualmente/Exclusivamente mutados  ####
#==============================================================================#
rop.si <- somaticInteractions(maf = maf.Tab, genes = genes.fdr$Hugo_Symbol, 
                    leftMar = 5, topMar = 5, 
                    pvalue = c(0.05, 0.05),fontSize = 0.6 )  
dev.off()

#==============================================================================#
# 9.10 Mutational Signatures      ####
#==============================================================================#
#Requires BSgenome object
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
rop.tnm = trinucleotideMatrix(maf = maf.Tab, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")

# 9.10.2 Differences between APOBEC enriched and non-enriched samples
#==============================================================================#
plotApobecDiff(tnm = rop.tnm, maf = maf.Tab, pVal = 0.2)

# 9.10.3 Signature analysis
#==============================================================================#
pdf(file = paste0(output,"Figures/signatures.pdf"), width = 12, height = 6)
#rop.sign = estimateSignatures(mat = rop.tnm, nTry = 8)
plotCophenetic(res = rop.sign)
rop.sig = extractSignatures(mat = rop.tnm, n = 5)
#View(rop.sig$contributions)

#Compate against original 30 signatures 
rop.og30.cosm = compareSignatures(nmfRes = rop.sig, sig_db = "legacy")
#Compate against updated version3 60 signatures 
rop.v3.cosm = compareSignatures(nmfRes = rop.sig, sig_db = "SBS")

pheatmap(mat = rop.og30.cosm$cosine_similarities, 
                   cluster_rows = FALSE, 
                   main = "cosine similarity against validated signatures")
plotSignatures(nmfRes = rop.sig, font_size = 1.5, title_size = 1.5, sig_db = "SBS")

dev.off()


library("barplot3d")
#Visualize first signature
sig1 = rop.sig$signatures[,1]
barplot3d::legoplot3d(contextdata = sig2, labels = T, scalexy = 0.1, 
                      sixcolors = "sanger", alpha = 0.5)


#==============================================================================#
#Survival analysis based on grouping of DNMT3A mutation status
#mafSurvival(maf = maf.Tab, genes = 'APC', time = 'ultimo.seg', Status = 'Follow.up', isTCGA = F)





#==============================================================================#
# Fisher Test  - recurrentemente mutated genes ####
# number of patients with a specific gene mutated, associated with the outcomes
#==============================================================================#
respostas <- c("Response1")

sink(file = paste0(output,"/Fisher_Test_dnds.mutGenes.FDR10.txt"), append = FALSE)
dndsGenes <-fread(paste0(output,"dndscv_pvalue_0.05.tsv"), header = T, na.strings=c("NA"), fill=TRUE, check.names = FALSE)
dndsGenes <- subset(dndsGenes, qglobal_cv <= 0.1, select = "gene_name")
for (resp in respostas) {
  for (gene in dndsGenes$gene_name) {
    aux <-Samples_Data
    aux$mutGene <-"not_dNdSGene"
    mutados<- unique(subset(syn.Mutation.MAF,  Gene.refGene == gene, select = Tumor_Sample_Barcode))
    aux[aux$Tumor_Sample_Barcode %in% mutados$Tumor_Sample_Barcode, "mutGene"] <-"dNdSGene"
    aux1<- as.data.frame(unique(subset(aux, select = c("Tumor_Sample_Barcode", "mutGene", resp))))
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
        print("========== Chisq Test =======================================")
        print(chisq.test(dt3))
        print("========== Fisher Test ======================================")
        print(fisher.test(dt3))
      }
    }
  }
}

sink()



