### Fig2B  #####################################################################
##  heatmap of tsRNA coverage
library(pheatmap)
library(tidyverse)
library(reshape2)

sample_name<-c("L3_ts1","L3_ts2","L3_ts3",
               "P1_ts1","P1_ts2",
               "P2_ts1","P2_ts2","P2_ts3")
sample_name1<-c("temper1","temper2","temper3")
sample_name2 <- c("a1","a2","a3","a4","a5","a6","a7","a8") 
# a1-a3: three replicats of L3; a4-a5: two repliactes of P1; a7-a9: three replicats of L3


help <- read.table("T2h.txt")
# In this file, each line placed alphabetically ranked tRNA isoacceptor as following.
# Ala-AGC
# Ala-CGC
# Ala-TGC
# ……


# Build a blank matrix
for(ic in 1:length(sample_name)){
  # four rows：tRNA.name start.position  end.position  counts 
  assign(sample_name1[1],read.table(paste(sample_name[ic],".txt",sep="")))
  assign(sample_name1[2],data.frame(matrix(0,43,74)))
  colnames(temper1) <- c("amino","start","length","count")
  temper1$end <- temper1$start+temper1$length-1
  temper1[(temper1$end >74),5] <- 74
  temper1[(temper1$start >74),2] <- 74
  colnames(temper2) <- c(1:74)
  rownames(temper2) <- help[,1]
  temper1[,4] <- temper1[,4]/sum(temper1[,4])*1000000
  # The reads location information was translated into the coverage of each site
  for(id in 1:nrow(temper2)){
    cc <- subset(temper1,amino == row.names(temper2)[id])
    for(ie in 1:nrow(cc)){
      temper2[id,cc[ie,2]:cc[ie,5]] <- temper2[id,cc[ie,2]:cc[ie,5]]+cc[ie,4]
    }
    temper2[id,] <- temper2[id,]/sum(cc[,4])*100
  }
  assign(sample_name2[ic],temper2) 
}


b1 <- a1
b2 <- a4
b3 <- a7
for(x in 1:43){for(y in 1:74){b2[x,y] <- (a4[x,y]+a5[x,y]+a6[x,y])/3}}
for(x in 1:43){for(y in 1:74){b1[x,y] <- (a1[x,y]+a2[x,y]+a3[x,y])/3}}
for(x in 1:43){for(y in 1:74){b3[x,y] <- (a7[x,y]+a8[x,y]+a9[x,y])/3}}
b4 <- rbind(b1,b2)
b4 <- rbind(b4,b3)
b5 <- b4
rn <- colnames(temper2)
rn[2:33] <- ""
rn[35:73] <- ""

pheatmap(b5,cluster_cols = F,cluster_rows = F,gaps_row = c(43,86),
         show_rownames=F,labels_col=rn,
         color = colorRampPalette(c("blue","#6666ff","white", "red"))(1000),border=FALSE)
