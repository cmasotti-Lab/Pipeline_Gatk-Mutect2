
library(data.table)
library(stringr)

wd <- "D:/Github-Projects/Pipeline_Gatk-Mutect2/Result_Mutect2.GATK4.6.2014-06-14/"
setwd(wd)
output <- "output-Step1/"



all<-fread(paste0(output,"all_variants.72sample.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)

#==============================================================================#
#             QUANTIFICANDO VARIANTES ####
#==============================================================================#
# 
## Quantificar variantes germinativas gnomAD and Abraom 
#==============================================================================#
nrow(all[all$gnomad41_exome_AF_grpmax > 0.005])  # 10165
nrow(all[all$gnomad41_exome_AF > 0.005])         # 7293
# Quantify mutation in HLA genes
#==============================================================================#
genesHLA <- fread("../input_Somatic-Filter/genesHLA.txt", col.names = "gene", header = F) 
all.hla<- subset(all, (Gene.refGene %in% genesHLA$gene))
nrow(unique(all.hla$index))    #95 mut em HLA
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
nrow(unique(all[shared2exclude,"index"]))  #10579 

# Quantificando mutações compartilhadas >=2, com menos de 4 registros no COSMIC?
#==============================================================================#
shared2exclude <- which((all$MUTATION_shared >= 2)  & (all$COSMIC_OCCURENCE <= 3))
nrow(unique(all[shared2exclude,"index"]))  # 2207 

nrow(unique(all[,"index"]))             # 79558
nrow(unique(all[,"SAMPLE"]))              # 72

dir.create(paste0(output,"Figures"))

aux1<-as.data.frame(table(all$COV), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
pdf(file = paste0(output,"Figures/Cov.all_variants.pdf"), width = 9, height = 5)

p<- ggplot(data=head(aux1, 500), aes(x=Var1, y=Freq))+
  geom_col() + 
  geom_vline(xintercept= 30, linetype="dashed", color = "red") +
  xlab("Deapth of Coverage") +theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
LEGEND<- c("LEGEND: The figure display the Coverage distribution range < 500x, from all variants in 72 samples.
           No one filter was applied. The redline indicates threashold at 20x.")
gridExtra::grid.arrange(p, bottom=LEGEND)
dev.off()

# Quantificando Distribuição  VAF 
#==============================================================================#
aux1<-as.data.frame(table(all$VAF), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
p<-ggplot(data=aux1, aes(x=Var1, y=Freq))+
  geom_col(width = .005) + 
  geom_vline(xintercept= 0.2, linetype="dashed", color = "red")+
  xlab("VAF") + theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
gridExtra::grid.arrange(p, bottom="Fig2: VAF distribution of variants in 72 samples")



#==============================================================================#
#    QUEM SÃO AS MUTAÇÕES COM VAF < 0.2 ?       ####

#==============================================================================#
## Quantifica MUTAÇOES not-rare gnomAD
#==============================================================================#
nrow(unique(all[all$gnomad41_exome_AF_grpmax >= 0.005, "index"]))
## Quantifica MUTAÇOES not-rare no Abraom
#==============================================================================#
nrow(unique(all[all$abraom_freq >= 0.005, "index"]))
## FILTRO: REMOVE mutações not-rare no gnomAD ou Abraom ####
#==============================================================================#
all.exclVAF <- all[c(not_rare2exclude,not_rareAbraom2exclude),]
nrow(unique(all.exclVAF$index))  # 5595

# Quantificando regiões das mutações
#==============================================================================#
aux<-(unique(all.exclVAF[,c("index", "Func.refGene")]))  
as.data.frame(table(aux$Func.refGene))
# Quantificando tipos de mutações
#==============================================================================#
aux<-(unique(all.exclVAF[,c("index", "ExonicFunc.refGene")]))  
as.data.frame(table(aux$ExonicFunc.refGene))
# Quantificando mutações compartilhadas >=2
#==============================================================================#
shared2exclude <- which((all.exclVAF$MUTATION_shared >= 2) )
nrow(unique(all.exclVAF[shared2exclude,"index"]))  #1221
# Quantificando mutações compartilhadas >=2, com menos de 4 registros no COSMIC?
#==============================================================================#
shared2exclude <- which((all.exclVAF$MUTATION_shared >= 2)  & (all.exclVAF$COSMIC_OCCURENCE <= 3))
nrow(unique(all.exclVAF[shared2exclude,"index"]))  # 680 

nrow(unique(all.exclVAF[,"index"]))             # 5595
nrow(unique(all.exclVAF$SAMPLE))              # 72

aux<-(unique(all.exclVAF[,c("index", "COV")]))  
summary(aux$COV)

aux1<-as.data.frame(table(all.exclVAF$COV), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
p<-ggplot(data=head(aux1, 1000), aes(x=Var1, y=Freq))+
  geom_col() + 
  geom_vline(xintercept= 10, linetype="dashed", color = "red")+theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
gridExtra::grid.arrange(p, bottom="COV distribution from 5.595 subclonal variants present in 72 samples")




#==============================================================================#
# SUMMARY SOMATIC VARIANTS
#==============================================================================#


all.x1<-as.data.frame(fread(paste0(output,"somatic_variants.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))
length(unique(all.x1$index)) # 73989 variantes somaticas

aux1<-as.data.frame(table(all.x1$VAF), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
p<-ggplot(data=aux1, aes(x=Var1, y=Freq))+
  geom_col(width = .005) + 
  geom_vline(xintercept= 0.2, linetype="dashed", color = "red")+ theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
gridExtra::grid.arrange(p, bottom="VAF distribution, 73.963 Mutations from 72 samples")




aux1<-as.data.frame(table(all.xg$COV), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
p<-ggplot(data=head(aux1, 500), aes(x=Var1, y=Freq))+
  geom_col() + 
  geom_vline(xintercept= 10, linetype="dashed", color = "red")+theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
gridExtra::grid.arrange(p, bottom="COV distribution, 64.956 Mutations from 72 samples")

aux1<-as.data.frame(table(all.xg$VAF), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
p<-ggplot(data=aux1, aes(x=Var1, y=Freq))+
  geom_col(width = .005) + 
  geom_vline(xintercept= 0.2, linetype="dashed", color = "red")+ theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
gridExtra::grid.arrange(p, bottom="VAF distribution, 64.956 Mutations from 72 samples")




somatic_Var.Filt_SC<-as.data.frame(fread(paste0(output,"somatic_variants.FiltSharedCosmic.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))

#==============================================================================#
#  QUANTIFICANDO MUTAÇÕES DO ALL MUTECT ANTES DOS FILTROS VAF E COV ####
#==============================================================================#
length(unique(somatic_Var.Filt_SC[,"index"]))             # 64978
length(unique(somatic_Var.Filt_SC$SAMPLE))              # 72

#  Quantify mutation in HLA genes
#==============================================================================#
genesHLA <- fread("../input_Somatic-Filter/genesHLA.txt", col.names = "gene", header = F) 
all.hla<- subset(somatic_Var.Filt_SC, (Gene.refGene %in% genesHLA$gene))
length(unique(all.hla$index))    # 55

# Quantificando regiões das mutações
#==============================================================================#
aux<-(unique(somatic_Var.Filt_SC[,c("index", "Func.refGene")]))  
as.data.frame(table(aux$Func.refGene))

# Quantificando tipos de mutações
#==============================================================================#
aux<-(unique(somatic_Var.Filt_SC[,c("index", "ExonicFunc.refGene")]))  
as.data.frame(table(aux$ExonicFunc.refGene))

# Quantificando Distribuição da VAF 
#==============================================================================#
aux<-(unique(somatic_Var.Filt_SC[,c("index", "VAF")]))  
summary(aux$VAF)
aux1<-as.data.frame(table(somatic_Var.Filt_SC$VAF), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
p<-ggplot(data=aux1, aes(x=Var1, y=Freq))+
  geom_col(width = .005) + 
  geom_vline(xintercept= 0.2, linetype="dashed", color = "red")+ theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
gridExtra::grid.arrange(p, bottom="VAF distribution, from 72 samples")
# Quantificando Distribuição da COV 
#==============================================================================#
aux<-(unique(somatic_Var.Filt_SC[,c("index", "COV")]))  
summary(aux$COV)
##PLOT distribuição da COV
#==============================================================================#
aux1<-as.data.frame(table(somatic_Var.Filt_SC$COV), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
p<-ggplot(data=head(aux1, 500), aes(x=Var1, y=Freq))+
  geom_col() + 
  geom_vline(xintercept= 10, linetype="dashed", color = "red")+theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
gridExtra::grid.arrange(p, bottom="COV distribution, from 85 samples")





#==============================================================================#
#    QUEM SÃO AS 34.464 MUTAÇÕES COM VAF < 0.1 ?       ####
#==============================================================================#
somatic_Var<-as.data.frame(fread(paste0(output,"somatic_variants.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))

VAF2exclude <- which(somatic_Var$VAF<0.10) 
length(unique(somatic_Var[VAF2exclude,"index"]))  #  34551
all.mut.VAF.excl <- somatic_Var[VAF2exclude,]
length(unique(all.mut.VAF.excl$index))  #19583

length(unique(all.mut.VAF.excl[,"index"]))             # 19583
length(unique(all.mut.VAF.excl$SAMPLE))              # 72

# Quantificando regiões das mutações
#==============================================================================#
aux<-(unique(all.mut.VAF.excl[,c("index", "Func.refGene")]))  
as.data.frame(table(aux$Func.refGene))
# Quantificando tipos de mutações
#==============================================================================#
aux<-(unique(all.mut.VAF.excl[,c("index", "ExonicFunc.refGene")]))  
as.data.frame(table(aux$ExonicFunc.refGene))
# Quantificando mutações compartilhadas >=2
#==============================================================================#
shared2exclude <- which((all.mut.VAF.excl$MUTATION_shared >= 2) )
length(unique(all.mut.VAF.excl[shared2exclude,"index"]))  #189
# Quantificando mutações compartilhadas >=2, com menos de 4 registros no COSMIC?
#==============================================================================#
shared2exclude <- which((all.mut.VAF.excl$MUTATION_shared >= 2)  & (all.mut.VAF.excl$COSMIC_OCCURENCE <= 3))
length(unique(all.mut.VAF.excl[shared2exclude,"index"]))  # 127 

length(unique(all.mut.VAF.excl[,"index"]))             # 19583
length(unique(all.mut.VAF.excl$SAMPLE))              # 72

aux1<-as.data.frame(table(all.mut.VAF.excl$VAF), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
p<-ggplot(data=aux1, aes(x=Var1, y=Freq))+
  geom_col(width = .005) + 
  geom_vline(xintercept= 0.10, linetype="dashed", color = "red")+ theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
gridExtra::grid.arrange(p, bottom="VAF distribution, from 72 samples")

aux<-(unique(all.mut.VAF.excl[,c("index", "COV")]))  
summary(aux$COV)

aux1<-as.data.frame(table(all.mut.VAF.excl$COV), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
p<-ggplot(data=head(aux1, 500), aes(x=Var1, y=Freq))+
  geom_col() + 
  geom_vline(xintercept= 10, linetype="dashed", color = "red")+theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
gridExtra::grid.arrange(p, bottom="COV distribution, 3.990 Mutations from 72 samples")






#==============================================================================#
# SUMMARY SOMATIC MUTATIONS , SYNONIMOUS INCLUED  ####
#==============================================================================#

somatic_Var.CodRegion<-as.data.frame(fread(paste0("output-Step1/somatic_variants.CodRegion.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))
length(unique(somatic_Var.CodRegion[,"index"]))  
# Quantificando mutações pelas regiões dos genes
#==============================================================================#
aux<-(unique(somatic_Var.CodRegion[,c("index", "Func.refGene")]))  
as.data.frame(table(aux$Func.refGene))

# Quantificando tipos de alterações
#==============================================================================#
aux<-(unique(somatic_Var.CodRegion[,c("index", "ExonicFunc.refGene")]))  
as.data.frame(table(aux$ExonicFunc.refGene))

aux<-(unique(somatic_Var.CodRegion[,c("index", "ExonicFunc.refGene")]))  
as.data.frame(table(aux$ExonicFunc.refGene))

# Generate a bar plot displaying distinct ExonicFunc mutation types
aux<-as.data.frame(unique(somatic_Var.CodRegion[,c("index", "ExonicFunc.refGene")]))
aux <- as.data.frame(table(aux$ExonicFunc.refGene))
aux <- aux %>%
  arrange(desc(Freq)) %>%
  mutate(Var1 = factor(Var1, levels = unique(Var1)))
colnames(aux) <- c("ExonicFunc", "count")
aux <- aux %>%  mutate(proportion = count / sum(count) * 100)
# Criar o gráfico com contagem e proporção
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

# Quantificando tipos de alterações
#==============================================================================#
length(unique(somatic_Var.CodRegion$Gene.refGene)) 

# Quantificando tipos de alterações
#==============================================================================#
p<- ggplot(data=sampleTableCount.Clinical, aes(x= reorder(EXOME_ID, -N_mutations), y=N_mutations, fill=Response1)) +
  geom_bar(stat="identity")+theme_classic() +
  scale_y_sqrt()+
  scale_fill_manual(values = c("#08415c", "#cc2936"), breaks = c("nCRT-R", "nCRT-NR")) +
  geom_hline(yintercept= median(sampleTableCount.Clinical$N_mutations), linetype="dashed", color = "black")+
  annotate(geom="text", x=40, y=500, label=paste("median N_mutations= ",sprintf("%.2f",median(sampleTableCount.Clinical$N_mutations))), color="black")+
  xlab("SAMPLES")+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="bottom") 
gridExtra::grid.arrange(p, bottom="Muatation distribution from 72 samples. *synonimous exclude")

p<- ggplot(data=sampleTableCount.Clinical, aes(x= reorder(EXOME_ID, -MATH_score), y=MATH_score, fill=Response1)) +
  geom_bar(stat="identity")+theme_classic() +
  scale_fill_manual(values = c("#08415c", "#cc2936"), breaks = c("nCRT-R", "nCRT-NR")) +
  geom_hline(yintercept= median(sampleTableCount.Clinical$MATH_score), linetype="dashed", color = "black")+
  annotate(geom="text", x=50, y=55, label=paste("median MATH_score= ",sprintf("%.2f",median(sampleTableCount.Clinical$MATH_score))), color="black")+
  xlab("SAMPLES")+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="bottom") 
gridExtra::grid.arrange(p, bottom="MATH_score distribution from 72 samples. *synonimous exclude")


#==============================================================================#
# DISTRIBUIÇÃO DAS VAFs
#==============================================================================#
aux1 <- select(subset(somatic_Var.CodRegion), SAMPLE, VAF)
aux2 <- as.data.frame(table(aux1), stringsAsFactors = F)
aux2$VAF <- as.double(aux2$VAF)
aux2 <- merge(aux2, Clinical[,c("Sample","EXOME_ID","Response1","Response2")], by.x= "SAMPLE", by.y = "Sample", all.x = T)

aux2 <- aux2[aux2$Freq >0,]

# Função para criar gráficos
create_plots <- function(data, response, color, output) {
  subset_data <- subset(data, Response2 == response)
  
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
  p3 <- ggplot(subset_data, aes(x = VAF, y = "", fill = Response2)) +
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
create_plots(aux2, "nCRT-NR Not-Metastatic", "#ef767a", output)
create_plots(aux2, "nCRT-NR Metastatic", "#cc2956", output)


