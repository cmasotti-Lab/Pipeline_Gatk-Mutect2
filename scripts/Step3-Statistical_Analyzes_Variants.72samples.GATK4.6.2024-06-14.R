
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
somatic<-as.data.frame(fread(("output-Step1/somatic_mutation.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))
somatic_Samples<- fread(paste0("output-Step1/MATH_TMB-SAMPLES.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)
somatic_Var.CodRegion<-as.data.frame(fread(paste0("output-Step1/somatic_variants.CodRegion.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))
sampleTableCount.Clinical<-fread(paste0(output,"/sampleTableCount.Clinical.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)
sampleTableCount.Clinical2<-fread(paste0(output,"/sampleTableCount.Clinical2.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)
driver<-fread("../input_Somatic-Filter/driver-genes.txt", header = F, col.names = "gene")
MMR<-fread("../input_Somatic-Filter/MMR_genes.txt", header = F, col.names = "gene")

#==============================================================================#
# Carregando tabela de dados Clinicos ####
#==============================================================================#
Clinical <- xlsx::read.xlsx("D:/PROJETOS-HSL_BP/Dados Clinicos 2024 ALLversions WXS-TCRseq.xlsx", sheetName = "73samples", header = T)
Clinical <- subset(Clinical, select=-Paciente)
Clinical <- Clinical[Clinical$WXS== "Y",]
Clinical <- Clinical[rowSums(is.na(Clinical)) != ncol(Clinical),]

#==============================================================================#
# DNDS_CV ####
#==============================================================================#

library(dndscv)

#==============================================================================#
# TESTE PARA CORRIGIR A ANALISE DE dNdS_cv PELA REGIÃO DO EXOMA
#==============================================================================#

# path_genome_fasta = system.file("extdata", "chr3_segment.fa", package = "dndscv", mustWork = TRUE)
# 
# target_file<- fread("../input_Somatic-Filter/xgen-exome-research-panel-v2-targets-hg38.bed ")
# target_file$Gene_v <-target_file$V4
# target_file <- separate(data = target_file, col = V4, into = c("V4"), sep = "_")
# target_file <- target_file[,c(4,4,7,1,2,3,2,3,3-2,6)]
# colnames(target_file) <- c("gene.id", "gene.name", "cds.id", "chr", "chr.coding.start", "chr.coding.end", "cds.start", "cds.end", "length", "strand")
# target_file[target_file$strand == "+",]$strand <- "1"
# target_file[target_file$strand == "-",]$strand <- "-1"
# fwrite(target_file, "../input_Somatic-Filter/xgen-exome-research-panel-v2-targets-hg38.cds", quote = F, col.names = T, row.names = F)
# fwrite(as.data.frame(unique(target_file$gene.name)), "../input_Somatic-Filter/xgen-exome-research-panel-v2-targets-hg38.list", quote = F, col.names = F, row.names = F)
# 

# cd input_Somatic-Filter/
# cut -f 4 xgen-exome-research-panel-v2-targets-hg38.bed | cut -d '_' -f1| sort -u > xgen-exome-research-panel-v2-targets-hg38.list
# path_cds_table <- c("../input_Somatic-Filter/xgen-exome-research-panel-v2-targets-hg38.mart_export.txt")
# path_cds_table <- c("../input_Somatic-Filter/mart_export.txt")
# path_cds_table <- c("../input_Somatic-Filter/chr19-hg38.mart_export.txt")
# 
# path_genome_fasta = c("../input_Somatic-Filter/Homo_sapiens.GRCh38.dna.chromosome.19.fa")
# path_genome_fasta = c("../input_Somatic-Filter/Homo_sapiens.GRCh38.p14_GCA_000001405.29.fa")
# buildref(cdsfile=path_cds_table, genomefile=path_genome_fasta, outfile = "example_output_refcds.rda", excludechrs="MT", useids = T,)

#==============================================================================#
# 72 SAMPLES
#==============================================================================#
all.mut.flt.72<-fread(paste0(output,"/somatic_variants.CodRegion.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)
mutations2 <- all.mut.flt.72[,c("SAMPLE", "Chr", "Start", "Ref", "Alt")]
mutations2$Chr <- gsub("chr", "", mutations2$Chr)
dndsout2 = dndscv(mutations2, refdb="D:/PROJETOS-HSL_BP/resultados_Mutect2-2023/input_Somatic-Filter/RefCDS_human_GRCh38_GencodeV18_recommended.rda", cv=NULL,
                  # gene_list = unique(targert_file$Gene),
                  max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf )


print(dndsout2$globaldnds)
# In particular, very low estimates of θ (the overdispersion parameter), 
# particularly θ<1, may reflect problems with the suitability of the dNdScv model for the dataset.
print(dndsout2$nbreg$theta)

signif_genes_localmodel = as.vector(dndsout2$sel_loc$gene_name[dndsout2$sel_loc$qall_loc<0.01])
print(signif_genes_localmodel)
# View(dndsout2$sel_cv) # This is shown as an example but these results based on a few genes should not be trusted
sel_cv = dndsout2$sel_cv
signif_genes = sel_cv[sel_cv$pglobal_cv<0.05, c("gene_name", "pglobal_cv" ,"qglobal_cv","qallsubs_cv")]
nrow(subset(signif_genes,pglobal_cv < 0.05))
nrow(subset(signif_genes,qglobal_cv < 0.1))
nrow(subset(signif_genes,qglobal_cv < 0.05))
nrow(subset(signif_genes,qallsubs_cv < 0.05)) # avaliando sem os Indels

fwrite(sel_cv[sel_cv$pglobal_cv< 0.05, ], paste0(output,"/dndscv_pvalue_0.05.tsv"), quote = F, sep="\t")
fwrite(sel_cv[sel_cv$qglobal_cv< 0.1, ], paste0(output,"/dndscv_FDR_0.1.tsv"), quote = F, sep="\t")
#==============================================================================#
#                ONCOPLOT  72 samples ####
#==============================================================================#
signif_genes <- fread(file=paste0(output,"/dndscv_FDR_0.1.tsv"))

somatic.VAF<-fread(paste0(output,"somatic_variants.CodRegion.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)
somatic.VAF<- subset(somatic.VAF, select=c("Chr", "Start",  "End","Ref", "Alt", "Gene.refGene","EXOME_ID","ExonicFunc.refGene"))

somatic.VAF<- somatic.VAF %>% 
  mutate(Gene.refGene = strsplit(as.character(Gene.refGene), ";")) %>% 
  unnest(Gene.refGene)

# DESCOMENTE ABAIXO PARA PLOTAR SEM MUTAÇOES SINONIMAS
# 

somatic.VAF$Reference_Allele <- somatic.VAF$Ref
colnames(somatic.VAF)[names(somatic.VAF)=="Ref"]<-"Tumor_Seq_Allele1"
colnames(somatic.VAF)[names(somatic.VAF)=="Alt"]<-"Tumor_Seq_Allele2"
colnames(somatic.VAF)[names(somatic.VAF)=="Gene.refGene"]<-"Hugo_Symbol"
colnames(somatic.VAF)[names(somatic.VAF)=="Chr"]<-"Chromosome"
colnames(somatic.VAF)[names(somatic.VAF)=="Start"]<-"Start_Position"
colnames(somatic.VAF)[names(somatic.VAF)=="End"]<-"End_Position"
colnames(somatic.VAF)[names(somatic.VAF)=="ExonicFunc.refGene"]<-"Variant_Classification"
somatic.VAF$Variant_Type<-"ExonicFunc.refGene"
colnames(somatic.VAF)[names(somatic.VAF)=="EXOME_ID"]<-"Tumor_Sample_Barcode"

sampleTableCount.Clinical.VAF <- sampleTableCount.Clinical
colnames(sampleTableCount.Clinical.VAF)[names(sampleTableCount.Clinical.VAF)=="EXOME_ID"]<-"Tumor_Sample_Barcode"

aux<-unique(subset(somatic, class == "MMR",select = c("SAMPLE","class") ))
sampleTableCount.Clinical.VAF$MMR = "pMMR"
sampleTableCount.Clinical.VAF$MMR[sampleTableCount.Clinical.VAF$SAMPLE %in% aux$SAMPLE] <-'MMR'

MAF = read.maf(maf = somatic.VAF,
               clinicalData = sampleTableCount.Clinical.VAF,
               vc_nonSyn = unique(somatic.VAF$Variant_Classification),
               verbose = F)
unique(somatic.VAF$Variant_Classification)


pdf(file = paste0(output,"Figures/oncoplot.pdf"), width = 12, height = 6)

oncoplot(maf = MAF,
         genes = subset(signif_genes, qglobal_cv <= 0.05, select=gene_name)$gene_name,
         sortByAnnotation = T,
         clinicalFeatures = c("Response1","Response2","Response3"), 
         showTumorSampleBarcodes = T, SampleNamefontSize = 0.8,
         removeNonMutated = F,  drawColBar = T,
         fontSize = 0.6, annoBorderCol = "white",
         annotationFontSize = 1, 
         legend_height =2, 
         legendFontSize = 1)

tmb(MAF, captureSize = 34.12676, logScale = TRUE)
dev.off()

#==============================================================================#
#  pheatmap - GENES with MMR  
#==============================================================================#
mut.genes=unique(select(somatic, c("Gene.refGene","EXOME_ID", "class")))
aux<- as.data.frame(unique(subset(mut.genes, select= EXOME_ID)))

mut.genes.heatmap= as.data.frame(table(unique(select(somatic, c("Gene.refGene","EXOME_ID")))))
mut.genes.heatmap<- mut.genes.heatmap[mut.genes.heatmap$Gene.refGene %in% MMR$gene,]
mut.genes.heatmap<- mut.genes.heatmap[mut.genes.heatmap$EXOME_ID %in% aux$EXOME_ID,]
df<- as.matrix(reshape2::dcast(mut.genes.heatmap, Gene.refGene ~ EXOME_ID,  value.var = "Freq"))
rownames(df) <- df[,1]
df<-type.convert(df[,-1], as.is = T)

## Juntando tabelas TMB_MATH_#MUT +  dados_clinicos 
my_sample_col <- merge(sampleTableCount, Clinical, by.x = "SAMPLE", by.y = "Sample")
my_sample_col <- data.frame(subset(Clinical, select = c("EXOME_ID","Response1","Response2","Response3")))

row.names(my_sample_col) <- my_sample_col[,1]
my_sample_col<-subset(my_sample_col, select=-1)

my_sample_col<-my_sample_col[order(my_sample_col$Response1 ),]
df_ordered <- df[,rownames(my_sample_col)]

pdf(file = paste0(output,"Figures/heatmap_dMMR.pdf"), width = 15, height = 8)

pheatmap((df_ordered), color=colorRampPalette(c("white", "darkred"))(100),
         annotation_col = my_sample_col,
         border_color = "grey60",
         cluster_rows  = F,cluster_cols  = F, fontsize = 8, 
         legend = F, main = "heatmap dMMR genes" )
dev.off()


#==============================================================================#
#        MUTALISK - ASSINATURAS MUTACIONAIS      ####
#==============================================================================#
#RODAR MUTALISK E BAIXAR PARA DENTRO DO DIRETORIO

# comando SHELL para cria a tabela
# grep -e '**THIS_SIGs\|**THIS_SIG_CONTRIBUTIONS\|**Cosine_Similarity' Mutalisk/*.txt | sed -E 's/_user.select_[0-9]_report.txt\:\*\*/\t/' | sed 's/Mutalisk\///' | sed 's/Cosine_Similarity\t/Cosine_Similarity/'| sed 's/\t/ /' |cut -d " " -f 1,2 > table_signatures.1.tsv
# grep -e '**THIS_SIGs\|**THIS_SIG_CONTRIBUTIONS\|**Cosine_Similarity' Mutalisk/*.txt | sed -E 's/_user.select_[0-9]_report.txt\:\*\*/\t/' | sed 's/Mutalisk\///' | sed 's/Cosine_Similarity\t/Cosine_Similarity/'| sed 's/\t/ /' |cut -d " " -f 3- | sed 's/ /;/g' > table_signatures.2.tsv

# grep -e '**THIS_SIGs\|**THIS_SIG_CONTRIBUTIONS\|**Cosine_Similarity' Mutalisk_synonymous/*.txt | sed -E 's/_user.selection_[0-9]_report.txt\:\*\*/\t/' | sed 's/Mutalisk_synonymous\///' | sed 's/Cosine_Similarity\t/Cosine_Similarity/'| sed 's/\t/ /' |cut -d " " -f 1,2 > table_signatures.1.synonymous.tsv
# grep -e '**THIS_SIGs\|**THIS_SIG_CONTRIBUTIONS\|**Cosine_Similarity' Mutalisk_synonymous/*.txt | sed -E 's/_user.selection_[0-9]_report.txt\:\*\*/\t/' | sed 's/Mutalisk_synonymous\///' | sed 's/Cosine_Similarity\t/Cosine_Similarity/'| sed 's/\t/ /' |cut -d " " -f 3- | sed 's/ /;/g' > table_signatures.2.synonymous.tsv


tab1<-fread(paste0("table_signatures.1.tsv"), header = F)
tab2<-fread(paste0("table_signatures.2.tsv"), sep="\t", header = F)

assinaturas<- cbind(tab1, tab2)
colnames(assinaturas) <- c("SAMPLE", "mutalisk", "value")
assinaturas<-dcast(assinaturas, SAMPLE ~ mutalisk)
assinaturas$SAMPLE <- gsub(".recode", "", assinaturas$SAMPLE)

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

final_Sig$MMR <-"pMMR"
mutados<- unique(subset(somatic, class == "MMR", select = SAMPLE))
final_Sig[final_Sig$SAMPLE %in% mutados$SAMPLE,"MMR"] <-"MMR"
# final_Sig<-merge(final_Sig, subset(somatic_Samples.resp, select = c("SAMPLE","Response1","Response3")), by="sample", all.y = T)

final_Sig<-merge(final_Sig, subset(sampleTableCount.Clinical, select = c("SAMPLE","Response1","Response3")), by="SAMPLE", all.y = T)

final_Sig$perc <-as.numeric(final_Sig$perc)
final_Sig$cosmic_sig <- "sig_others"
final_Sig[final_Sig$sig_id %in% c("6","15","20","26"),"cosmic_sig"] <-"MMR_sig"

fwrite(final_Sig, "table_signatures.tsv", sep = "\t", quote = F, col.names = T, row.names = F)

mutsig <-fread("table_signatures.tsv", header = T)
mutsig <- merge(Clinical[,c(1,2)], mutsig, by.x= "Sample",by.y = "SAMPLE")
#mutsig<-subset(mutsig, MMR == "MMR")
mutsig$sig_id<-as.character(mutsig$sig_id)

length(unique(mutsig$EXOME_ID))
length(unique(mutsig$Sample))
mutsig <- subset(mutsig, EXOME_ID %in% rownames(my_sample_col))
mutsig$EXOME_ID <- factor(mutsig$EXOME_ID, levels=rownames(my_sample_col))
my_sig_id<-sort(as.numeric(unique(mutsig$sig_id)))

pdf(file = paste0(output,"Figures/mutalisk_cosmicv2.pdf"), width = 12, height = 8)

p<-ggplot(data=mutsig, aes(x=EXOME_ID, y=perc, fill=sig_id)) +
  geom_bar(position="stack", stat="identity") + theme_classic() +xlab("")+
  scale_fill_manual(name="Mutational Signatures Cosmic_v2",
                    #breaks = my_sig_id,
                    breaks = c("10","6","15","20", "3","1","2","4","5","7","8","9","12","13","14","16","17","18","21","22","24","28","11","30"),
                    values=c("#f3722c","forestgreen","darkgreen","yellowgreen", "#edd3f9", "#edE5f9", "#edE2f9","#edf2f4","#fefae0","#edede9","#d6ccc2","#e3d5ca","#e3d45a","#fec5bb","#fae1dd","#c9ada7","#ffc999","#ffcdb2","#ffd7ba","#fec89a","#d5bdaf","#cad2c5","#edede9","#d6ccc2","#e3d5ca","#e3d5ca" )
  )
p +theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=0.2), legend.position = "bottom")

dev.off()





library(ggplot2)
library(data.table)
library(dplyr)
library(scales)

# Função para processar e gerar gráficos
process_and_plot <- function(file_path, clinical_data, output, plot_title) {
  # Carregar dados
  data <- fread(file_path, header = TRUE)
  
  # Transformar e processar dados
  data <- data %>%
    melt(id.vars = c("Signature", "Proposed_Etiology"), variable.name = "SAMPLE", value.name = "perc") %>%
    mutate(SAMPLE = gsub("\\.", "-", SAMPLE)) %>%
    merge(subset(clinical_data, select = c("Sample", "EXOME_ID", "Response1", "Response2")), 
          by.x = "SAMPLE", by.y = "Sample") %>%
    filter(perc > 0)
  
  # Gerar cores para as assinaturas
  num_signatures <- length(unique(data$Signature))
  signature_colors <- hue_pal()(num_signatures)
  
  # Criar o gráfico e salvar em PDF
 # pdf(file = paste0(output,"/Figures/", plot_title, ".pdf"), width = 12, height = 8)
  ggplot(data, aes(x = EXOME_ID, y = perc, fill = Signature)) +
    geom_bar(position = "stack", stat = "identity") +
    theme_classic() +
    xlab("") +
    scale_fill_manual(name = plot_title, values = signature_colors) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 0.2), 
          legend.position = "bottom") +
    labs(title = plot_title)
#  dev.off()
}

pdf(file = paste0(output,"/Figures/SigProfiles.pdf"), width = 12, height = 8)
# Processar e gerar gráficos para cada conjunto de dados
process_and_plot("SigProfiler_SBS/Signature_contributions.csv", Clinical, output, "SigProfile SBS")
process_and_plot("SigProfiler_DBS/Signature_contributions.csv", Clinical, output, "SigProfile DBS")
process_and_plot("SigProfiler_ID/Signature_contributions.csv", Clinical, output, "SigProfile ID")
dev.off()



#==============================================================================#
# Fisher Test  - recurrentemente mutated genes ####
# number of patients with a specific gene mutated, associated with the outcomes
#==============================================================================#
respostas <- c("Response1", "Response2", "Response3")

sink(file = paste0(output,"/Fisher_Test_dnds.mutGenes.FDR10.txt"), append = FALSE)
dndsGenes <-fread(paste0(output,"dndscv_pvalue_0.05.tsv"), header = T, na.strings=c("NA"), fill=TRUE, check.names = FALSE)
dndsGenes <- subset(dndsGenes, qglobal_cv <= 0.1, select = "gene_name")
for (resp in respostas) {
  for (gene in dndsGenes$gene_name) {
    aux <-sampleTableCount.Clinical
    aux$mutGene <-"not_dNdSGene"
    mutados<- unique(subset(somatic,  Gene.refGene == gene, select = EXOME_ID))
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
        print("========== Chisq Test =======================================")
        print(chisq.test(dt3))
        print("========== Fisher Test ======================================")
        print(fisher.test(dt3))
      }
    }
  }
}

sink()
