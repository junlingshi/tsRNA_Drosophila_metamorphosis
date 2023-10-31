###The five columns in P2_L3_g are L3_1, L3_2 L3_3, P2_1 and P2_2;
library(DEseq2)
rownames(P2_L3_g)<-P2_L3_g$Geneid
P2_L3_g<-P2_L3_g[,-1]
colData<-data.frame(row.names=colnames(P2_L3_g),condition)
dds<-DESeqDataSetFromMatrix(P2_L3_g,colData,design=~condition)
dds<-DESeq(dds)
P2_L3_g<-results(dds)
P2_L3_g<-as.data.frame(P2_L3_g)
P2_L3_g<-subset(P2_L3_g,baseMean>70)
P2_L3_g$Geneid<-rownames(P2_L3_g)
P2_L3_g<-merge(P2_L3_g,target_15,by="Geneid",all.x=T)
P2_L3_g<-unique(P2_L3_g)
ggplot(P2_L3_g,aes(type,log2FoldChange,fill=type))+geom_boxplot()+theme_bw()+
  scale_x_discrete(limits=c("No sites","others","meta","trans"))+
  geom_signif(comparisons = list(c("others","No sites"),c("meta","No sites"),c("trans","No sites")),map_signif_level = T,textsize = 4,test = wilcox.test,step_increase = 0.2)+theme(text=element_text(size=14,  family="serif"))
ggplot(P2_L3_g,aes(tar,log2FoldChange,fill=tar))+geom_boxplot()+theme_bw()+
  scale_x_discrete(limits=c("No sites","targets"))+
  geom_signif(comparisons = list(c("targets","No sites")),map_signif_level = T,textsize = 4,test = ks.test,step_increase = 0.2)+theme(text=element_text(size=14,  family="serif"))+coord_cartesian(xlim=c(-2,2))
