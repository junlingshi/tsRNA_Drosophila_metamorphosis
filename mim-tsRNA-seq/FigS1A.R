### FigS1A ####################################################################
aa <- read.table("tsRNA.txt",header = T)
aa[,2] <- as.numeric(aa[,2])
aa[,3] <- as.numeric(aa[,3])
aa[,4] <- as.numeric(aa[,4])
aa <- subset(aa,abs((end-start+1)-25.5) <= 3.5)
colnames(aa) <- c("tRNA","start","end","count")
length <- read.table("length.txt")
colnames(length) <- c("tRNA","length")
aa <- merge(aa,length,by="tRNA")
aa$length <- aa$length+3

aa$start <- aa$start/aa$length*100
aa$end <- aa$end/aa$length*100
aa$start <- round((aa$start/4)-0.000001+0.5)
aa$end <- round((aa$end/4)-0.000001+0.5)
aa <- aa[,-5]
aa <- subset(aa,tRNA != "SeC")


tRNA <- as.data.frame(matrix(ncol=59,nrow=length(unique(aa$tRNA))))
a <- c("tRNA")
a[2:59] <- c(1:58)
colnames(tRNA) <- a
b <- unique(aa$tRNA)
tRNA[,1] <- b
tRNA[,2:59] <- 0
i <- 1
for(x in 1:nrow(aa)){
  if (aa[x,1] != b[i]) {
    i <- i+1
    tRNA[i,(aa[x,2]+1):(aa[x,3]+1)] <- tRNA[i,(aa[x,2]+1):(aa[x,3]+1)]+aa[x,4]
  }
  else {
    tRNA[i,(aa[x,2]+1):(aa[x,3]+1)] <- tRNA[i,(aa[x,2]+1):(aa[x,3]+1)]+aa[x,4]
  }
}


tRNA_1 <- tRNA[,1:26]
tRNA_1 <- subset(tRNA_1,tRNA != "SeC-TCA")
d <- 0
for(x in 2:26){d[x-1] <- sum(tRNA_1[,x])}
tRNA_2 <- tRNA_1[order(tRNA_1$`9`),]
for(x in 2:26){tRNA_2[,x] <- tRNA_2[,x]/sum(tRNA_1[,9])}
c <- tRNA_2$tRNA

tRNA_2 <- melt(tRNA_2,id="tRNA")
tRNA_2$variable <- factor(tRNA_2$variable)
#tRNA_2$tRNA <- factor(tRNA_2$tRNA,levels=c)
tRNA_2$tRNA <- factor(tRNA_2$tRNA,levels=rate)

colour <- c("#b8d8e8","#87b9d7","#5699c6","#75b0b6","#a9d4a8","#a8d98d",
            "#79c16c","#74b262","#bcb089","#f9a8a8","#f1797b","#e94c4d",
            "#e94c4d","#f07d63","#f9b881","#fdbc71","#fda549","#e3b2a1",
            "#cdb8d9","#a98dc3","#8863ae")

colour1 <- colorRampPalette(colors = colour)(44)
library(ggsci)
pic <- ggplot(tRNA_2,aes(x = variable, y=value, fill=tRNA, group= tRNA))
pic <- pic + 
  geom_bar(stat="identity") + geom_col(position = 'stack', width = 0.6)+ 
  labs(y = "Scaled normalized coverage")+
  theme_bw()+scale_fill_manual(values=colour1)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

show(pic)

graph2ppt(pic,file="pic/tRNA_cover_ratio_shortreads_37.pptx",height=5,width=9)
