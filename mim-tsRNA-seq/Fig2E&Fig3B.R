L3_class1<-read.table(file="L3-ts1_class.txt")
L3_class2<-read.table(file="L3-ts2_class.txt")
L3_class3<-read.table(file="L3-ts3_class.txt")
P1_class1<-read.table(file="P1-ts1_class.txt")
P1_class2<-read.table(file="P1-ts2_class.txt")
P2_class1<-read.table(file="P2-ts1_class.txt")
P2_class2<-read.table(file="P2-ts2_class.txt")
P2_class3<-read.table(file="P2-ts3_class.txt")
L3_class_RPKR<-merge(L3_class1,L3_class2,by="V1")
L3_class_RPKR<-merge(L3_class_RPKR,L3_class3,by="V1")
colnames(L3_class_RPKR)<-c("tRNA_type","L3_ts1","L3_ts2","L3_ts3")
L3_class_RPKR$stage<-c("L3")
P1_class_RPKR<-merge(P1_class1,P1_class2,by="V1")
P2_class_RPKR<-merge(P2_class1,P2_class2,by="V1")
P2_class_RPKR<-merge(P2_class_RPKR,P2_class3,by="V1")
colnames(P2_class_RPKR)<-c("tRNA_type","P2_ts1","P2_ts2","P2_ts3")
colnames(P1_class_RPKR)<-c("tRNA_type","P1_ts1","P1_ts2")
P1_class_RPKR$stage<-c("P1")
P2_class_RPKR$stage<-c("P2")
L3_class_RPKR$L3_ts1<-L3_class_RPKR$L3_ts1/570012*1000
L3_class_RPKR$L3_ts2<-L3_class_RPKR$L3_ts2/153358*1000
L3_class_RPKR$L3_ts3<-L3_class_RPKR$L3_ts3/308084*1000
P1_class_RPKR$P1_ts1<-P1_class_RPKR$P1_ts1/106290*1000
P1_class_RPKR$P1_ts3<-P1_class_RPKR$P1_ts2/233521*1000
P2_class_RPKR$P2_ts1<-P2_class_RPKR$P2_ts1/684997*1000
P2_class_RPKR$P2_ts2<-P2_class_RPKR$P2_ts2/130698*1000
P2_class_RPKR$P2_ts3<-P2_class_RPKR$P2_ts3/434355*1000
P2_class_RPKR$RPKR<-(P2_class_RPKR$P2_ts3+P2_class_RPKR$P2_ts2+P2_class_RPKR$P2_ts1)/3
L3_class_RPKR$RPKR<-(L3_class_RPKR$L3_ts3+L3_class_RPKR$L3_ts2+L3_class_RPKR$L3_ts1)/3
P1_class_RPKR$RPKR<-(P1_class_RPKR$P1_ts2+P1_class_RPKR$P1_ts1)/2
L3_class_RPKR<-separate(L3_class_RPKR,tRNA_type,into = c("tRNA","type"),sep="_")
P1_class_RPKR<-separate(P1_class_RPKR,tRNA_type,into = c("tRNA","type"),sep="_")
P2_class_RPKR<-separate(P2_class_RPKR,tRNA_type,into = c("tRNA","type"),sep="_")
L3_RPKR<-L3_class_RPKR[,c("tRNA","type","RPKR","stage")]
P1_RPKR<-P1_class_RPKR[,c("tRNA","type","RPKR","stage")]
P2_RPKR<-P2_class_RPKR[,c("tRNA","type","RPKR","stage")]
RPKR<-rbind(L3_RPKR,P1_RPKR,P2_RPKR)
library(reshape2)
ff<-melt(RPKR)
ff<-ff[,-4];ff<-dcast(ff,tRNA~type+stage);rownames(ff)<-ff$tRNA
ee<-ff[,c("5-tRF_L3","5-tRF_P1","5-tRF_P2","5-tRH_L3","5-tRH_P1","5-tRH_P2","Inter-tRF_L3","Inter-tRF_P1","Inter-tRF_P2","3-tRH_L3","3-tRH_P1","3-tRH_P2","3-tRF_L3","3-tRF_P1","3-tRF_P2")]
library(pheatmap)
pheatmap(log(ee,2),cluster_cols = F,angle_col =90,fontsize_col = 10,gaps_col =c(3,6,9,12),main = "nuclearly encoded tRNA",color = c(colorRampPalette(colors = c("blue","white","red"))(100)))
