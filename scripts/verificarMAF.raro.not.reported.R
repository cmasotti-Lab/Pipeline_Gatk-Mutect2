

tab <- all.x1

table(tab$freq_gnomAD_lower01)
table(tab$freq_abraom_status)
#==============================================================================#
 raros <-subset(tab,  freq_gnomAD_lower01 =="rare" )
# raros1 <- which(tab$freq_gnomAD_lower01=="rare")
# raros2 <- which(tab$freq_abraom_status=="rare")
# raros <- tab[c(raros1,raros2),]

aux<-(unique(raros[,c("index", "MAF")]))  
summary(aux$MAF)
aux1<-as.data.frame(table(raros$MAF), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
p<-ggplot(data=aux1, aes(x=Var1, y=Freq))+
  geom_col(width = .005) + 
  geom_vline(xintercept= 0.2, linetype="dashed", color = "red")+ theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
gridExtra::grid.arrange(p, bottom="rare")
 
table(tab$freq_gnomAD_lower01)
table(tab$freq_abraom_status)

rareMAF5.1 <- subset(all.x1,  MAF == 1)
barplot(table(rareMAF5.1$Chr))

teste<-(rareMAF5.1[rareMAF5.1$SAMPLE %in% Clinical[Clinical$Gender == "M",]$Sample,])
barplot(table(teste$Chr))
#==============================================================================#
Nraros <-subset(tab,  freq_gnomAD_lower01 =="not-reported" )
# Nraros1 <- which(tab$freq_gnomAD_lower01 =="not-reported")
# Nraros2 <- which(tab$freq_abraom_status =="not-reported")
Nraros <- tab[c(Nraros1, Nraros2),]

aux<-(unique(Nraros[,c("index", "MAF")]))  
summary(aux$MAF)
aux1<-as.data.frame(table(Nraros$MAF), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
p<-ggplot(data=aux1, aes(x=Var1, y=Freq))+
  geom_col(width = .005) + 
  geom_vline(xintercept= 0.2, linetype="dashed", color = "red")+ theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
gridExtra::grid.arrange(p, bottom="not-reported")


#==============================================================================#
aux<-(unique(all.mut.maf.excluded[,c("index", "COV")]))  
summary(aux$COV)
aux1<-as.data.frame(table(all.mut.maf.excluded$COV), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
p<-ggplot(data=aux1, aes(x=Var1, y=Freq))+
 theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
gridExtra::grid.arrange(p, bottom="cov")


aux1<-as.data.frame(table(all.mut.maf.excluded$MAF), stringsAsFactors = F)
aux1$Var1<-as.numeric(aux1$Var1)
p<-ggplot(data=aux1, aes(x=Var1, y=Freq))+
  geom_col(width = .005) + 
  geom_vline(xintercept= 0.2, linetype="dashed", color = "red")+ theme_classic()+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")  
gridExtra::grid.arrange(p, bottom="MAF < 20")

boxplot(all.mut.maf.excluded$COV)
hist(all.mut.maf.excluded$COV)
