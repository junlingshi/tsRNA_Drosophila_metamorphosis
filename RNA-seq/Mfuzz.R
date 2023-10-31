tt<-subset(mfuzz_TE,L3_mRNA>5 & P1_mRNA >5 & P2_mRNA >5)
tt<-subset(tt,abs(log(P1_L3_FC,2))>1 | abs(log(P2_P1_FC,2))>1)
cc<-tt[,c("Geneid","L3_TE","P1_TE","P2_TE")]
row.names(cc)<-cc$Geneid
cc<-cc[,-1]
cc<-as.matrix(cc)
aa<-new('ExpressionSet',exprs=cc)
aa<-filter.NA(aa,thres=0.25)
aa<-fill.NA(aa,mode='mean')
aa<-filter.std(aa,min.std = 0)
aa<-filter.std(aa,min.std = 0)
aa<-standardise(aa)
c<-4
cl<-mfuzz(aa,c=c,m=m)
color.2 <- colorRampPalette(rev(c("#949494", "#a9a9a9", "#bdbebd")))(10000)
mfuzz.plot(aa,cl=cl,mfrow = c(1,4),colo = color.2)
gene_cluster <- cl$cluster
cl$size
gene_standard <- aa@assayData$exprs
gene_standard_cluster <- cbind(gene_standard[names(gene_cluster), ], gene_cluster)
gene_standard_cluster<-as.data.frame(gene_standard_cluster)
gene_standard_cluster$Geneid<-rownames(gene_standard_cluster)
gene_standard_cluster<-merge(gene_standard_cluster,target_e,by="Geneid")
meta_mRNA2<-subset(gene_standard_cluster,gene_cluster=="2" & type=="meta")
meta_mRNA2<-meta_mRNA2[,c("Geneid","L3","P1","P2")]
meta_mRNA2<-melt(meta_mRNA2)
ggplot(meta_mRNA2,aes(variable,value,group=Geneid))+geom_point()+geom_line()+theme_bw()
