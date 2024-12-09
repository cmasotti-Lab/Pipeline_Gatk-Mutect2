library(maftools)




wd <- "D:/Github-Projects/Pipeline_LAPC/"
setwd(wd)

# Create output directory with current date
#output <- paste0("OUTPUT-Somatic_Filter.",Sys.Date())
output <- paste0("OUTPUT-Somatic_Filter.2024-06-13")
dir.create(output)

rop_all <- fread(file="Final_GATK.4.6_annotated.txt", header=T, sep="\t", na.strings=c(".","NA"), fill=TRUE, check.names = FALSE)
rop.all.maftoos = annovarToMaf(annovar = "Final_GATK.4.6_annotated.txt", Center = 'CSI-NUS', refBuild = 'hg38', 
                           tsbCol = 'Tumor_Sample_Barcode', table = 'refGene', )
somatic_Var.CodRegion<-as.data.frame(fread(paste0(output,"somatic_variants.CodRegion.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))

rop.somatic = annovarToMaf(annovar = "output-Step1/somatic_variants.CodRegion.tsv", Center = 'CSI-NUS', refBuild = 'hg38', 
                     tsbCol = 'Tumor_Sample_Barcode', table = 'refGene', )

dim(somatic_Var.CodRegion)
dim(rop.somatic)


aux.colnames <- colnames(rop.somatic[,c(14:20)])
aux1<- subset(rop.somatic, select = aux.colnames)
aux2<- subset(somatic_Var.CodRegion, select = aux.colnames)
aux3<-setdiff(aux1, aux2)

aux <- subset(rop.somatic, select = c(8, 18,20))

rop.1 = annovarToMaf(annovar = "ROP-1-ExC85-xgenV2_S7.hg38_multianno.txt", Center = 'CSI-NUS', refBuild = 'hg38', 
                               tsbCol = 'Tumor_Sample_Barcode', table = 'refGene', )
rop.all = annovarToMaf(annovar = "Final_GATK.4.6_annotated.txt", 
                       Center = 'ROP', refBuild = 'hg38', 
                       tsbCol = 'Tumor_Sample_Barcode', table = 'refGene', )

var.annovar.maf$Tumor_Sample_Barcode = var.annovar.maf$sample_id
maf.Tab <- read.maf(maf = var.annovar.maf,
                    clinicalData = clinicalData.samples,
                    vc_nonSyn = unique(var.annovar.maf$Variant_Classification),
                    verbose = F)
tcga.ab.2972.het = inferHeterogeneity(maf = maf.Tab, 
                                      tsb = 'ROP-1-ExC85-xgenV2_S7', 
                                      vafCol = 'i_TumorVAF_WU',)

laml.mutsig.corrected = prepareMutSig(maf = maf.Tab)


Tab.toMATH <- somatic_Var.CodRegion

Tab.toMATH <- subset(Tab.toMATH, select=c("Chr", "Start", "End","Ref", "Alt", "Gene.refGene","EXOME_ID","ExonicFunc.refGene", "AAChange.refGene", "VAF"))
#Tab.toMATH <- subset(Tab.toMATH$, select=c("Chr", "Start", "End","Ref", "Alt", "Gene.refGene","EXOME_ID","ExonicFunc.refGene", "VAF"))

# DESCOMENTE ABAIXO PARA PLOTAR SEM MUTAÇOES SINONIMAS
Tab.toMATH$Reference_Allele <- Tab.toMATH$Ref
colnames(Tab.toMATH)[names(Tab.toMATH)=="Ref"]<-"Tumor_Seq_Allele1"
colnames(Tab.toMATH)[names(Tab.toMATH)=="Alt"]<-"Tumor_Seq_Allele2"
colnames(Tab.toMATH)[names(Tab.toMATH)=="Gene.refGene"]<-"Hugo_Symbol"
colnames(Tab.toMATH)[names(Tab.toMATH)=="Chr"]<-"Chromosome"
colnames(Tab.toMATH)[names(Tab.toMATH)=="Start"]<-"Start_Position"
colnames(Tab.toMATH)[names(Tab.toMATH)=="AAChange.refGene"]<-"AAChange"
colnames(Tab.toMATH)[names(Tab.toMATH)=="End"]<-"End_Position"
colnames(Tab.toMATH)[names(Tab.toMATH)=="ExonicFunc.refGene"]<-"Variant_Classification"
Tab.toMATH$Variant_Type<-"ExonicFunc.refGene"
colnames(Tab.toMATH)[names(Tab.toMATH)=="EXOME_ID"]<-"Tumor_Sample_Barcode"

clinicalData.samples <- as.data.frame(unique(Tab.toMATH$Tumor_Sample_Barcode))
colnames(clinicalData.samples) <- "Tumor_Sample_Barcode"

respondedores <- Clinical[Clinical$Response1 == "nCRT-R", "EXOME_ID"]
Tab.toMATH$Status = 0


primary.apl <- subset(Tab.toMATH, Tumor_Sample_Barcode %in% Clinical[Clinical$Response1 == "nCRT-R", "EXOME_ID"], )
relapse.apl <- subset(Tab.toMATH, Tumor_Sample_Barcode %in% Clinical[Clinical$Response1 == "nCRT-NR", "EXOME_ID"], )
  
primary.apl <- read.maf(maf = primary.apl,
                        clinicalData = clinicalData.samples,
                        vc_nonSyn = unique(Tab.toMATH$Variant_Classification),
                        verbose = F)
relapse.apl <- read.maf(maf = relapse.apl,
                        clinicalData = clinicalData.samples,
                        vc_nonSyn = unique(Tab.toMATH$Variant_Classification),
                        verbose = F)
pt.vs.rt <- mafCompare(m1 = primary.apl, m2 = relapse.apl, m1Name = 'Primary',
                       m2Name = 'Relapse', minMut = 5)
View(pt.vs.rt$results)


dnds_genes <- fread("output-Step1/dndscv_mutGENES.FDR_0.1.tsv")
dnds_genes <- dnds_genes[dnds_genes$qglobal_cv <= 0.1,]$gene_name
maf.Tab <- read.maf(maf = subset(Tab.toMATH,Hugo_Symbol %in% dnds_genes  ),
                    clinicalData = clinicalData.samples,
                    vc_nonSyn = unique(Tab.toMATH$Variant_Classification),
                    verbose = F)


# Calcular interações somáticas (co-ocorrência e exclusividade mútua)
# O parâmetro 'top' define o número de genes mais frequentemente mutados a serem considerados
pdf(file = paste0(output,"co-ocorrencia.pdf"), width = 9, height = 5)
somatic_interactions <- somaticInteractions(maf = maf.Tab, top = 20, pvalue = c(0.05, 0.1))
dev.off()

# Alternativamente, você pode usar a função oncodrive
oncodrive_results <- oncodrive(maf = maf.Tab, minMut = 5, pvalMethod = 'zscore')
View(oncodrive_results)

mut_freq <- mafSurvival(maf = maf_data, genes = NULL, time = "Time", Status = "Status", isTCGA = FALSE)



View(laml@data)


laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read.maf(maf = laml.maf)
tmb(maf = laml,  )
