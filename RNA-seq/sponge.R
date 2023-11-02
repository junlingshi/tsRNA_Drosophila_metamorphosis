setwd("D://data/Drosophila tsRNA/Seq/RNAseq/sponge/")
Gly<-read.table(file = "Gly_count.txt",header = TRUE)
P2<-read.table(file = "P2_Scr.txt")
Gly<-merge(P2,Gly,by="Geneid")
rownames(Gly)<-Gly$Geneid
Gly<-Gly[,-1]
colnames(Gly)<-c("P2-1","P2-2","Gly-1","Gly-2")
condition <- factor(c("control","control","treat","treat"), levels = c("control","treat"))
colData<-data.frame(row.names=colnames(Gly),condition)
library(DESeq2)
dds<-DESeqDataSetFromMatrix(Gly,colData,design = ~condition)
dds<-DESeq(dds)
Gly_d<-as.data.frame(results(dds))
Gly_d$Geneid<-rownames(Gly_d)

###The influence of tsRNA site number on target genes repression effect
Gly_d$Geneid<-rownames(Gly_d)
Gly_target<-read.table(file="Gly_target_0.05.txt")
colnames(Gly_target)<-c("Geneid","type")
c<-setdiff(Gly_d$Geneid,Gly_target$Geneid)
no_Gly<-Gly_d[which(Gly_d$Geneid %in%c),]
no_Gly$number<-c("No sites")
no_Gly$site_num<-c("No sites")
Gly_site_num<-read.table(file="Gly_site_num_s_0.05.txt")
colnames(Gly_site_num)<-c("Geneid","number")
Gly_site_num$number<-as.numeric(Gly_site_num$number)
Gly_s<-merge(Gly_d,Gly_site_num,by="Geneid")
Gly_s$site_num<-as.factor(ifelse(Gly_s$number>=1 & Gly_s$number<=2,'1-2 sites',ifelse(Gly_s$number>=3 & Gly_s$number<=4,'3-4 sites',ifelse(Gly_s$number==0,'No sites','>=5 sites'))))
Gly_s<-rbind(no_Gly,Gly_s)
Gly_s<-Gly_s[complete.cases(Gly_s),]
Gly_s<-subset(Gly_s,baseMean>50)
ggplot(Gly_s,aes(log2FoldChange))+stat_ecdf(aes(group=site_num,color=site_num))+theme_bw()+xlab("Gly")+coord_cartesian(xlim =c(-2,2))
ggplot(Gly_s,aes(site_num,log2FoldChange,fill=site_num))+geom_boxplot()+theme_bw()+xlab("Gly")+theme(text=element_text(size=14,  family="serif"))+scale_x_discrete(limits=c("No sites","1-2 sites","3-4 sites",">=5 sites"))+geom_signif(comparisons = list(c("1-2 sites","No sites"),c("3-4 sites","No sites"),c(">=5 sites","No sites")),map_signif_level = T,textsize = 4,test = ks.test,step_increase = 0.2)

### analysis of tsRNA target site position biase
Gly_3<-read.table(file="Gly_FBgn_3UTR.txt")
Gly_3$type<-c("3'UTR")
colnames(Gly_3)<-c("Geneid","type")
Gly_5<-read.table(file="Gly_FBgn_5UTR.txt")
Gly_5$type<-c("5'UTR")
colnames(Gly_5)<-c("Geneid","type")
Gly_C<-read.table(file="Gly_FBgn_CDS.txt")
Gly_C$type<-c("CDS")
colnames(Gly_C)<-c("Geneid","type")
Gly_target<-rbind(Gly_3,Gly_5,Gly_C)
Gly_d$Geneid<-rownames(Gly_d)
c<-setdiff(Gly_d$Geneid,Gly_target$Geneid)
no_Gly<-Gly_d[which(Gly_d$Geneid %in%c),]
no_Gly$type<-c("No sites")
Gly_t<-merge(Gly_d,Gly_tar_uni,by="Geneid")
Gly_f<-rbind(Gly_t,no_Gly)
Gly_f<-Gly_f[complete.cases(Gly_f),]
Gly_f<-subset(Gly_f,baseMean>50)
ggplot(Gly_f,aes(type,log2FoldChange,fill=type))+geom_boxplot()+theme_bw()+
  scale_x_discrete(limits=c("No sites","5'UTR","CDS","3'UTR"))+
  geom_signif(comparisons = list(c("5'UTR","No sites"),c("CDS","No sites"),c("3'UTR","No sites")),map_signif_level = T,textsize = 4,test = t.test,step_increase = 0.2)+
  xlab("Gly")+theme(text=element_text(size=14,  family="serif"))
